package ApproxSS

import (
	"fmt"
	"io"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// FastShamirApproxSS 基于MyShamirApproxSS改进的算法，使用加法全同态加密和门限解密
type FastShamirApproxSS struct {
	VanSS *VanillaShamirSS

	params4DoubleEncryption ckks.Parameters
	params4cmb              rlwe.Parameters
	f                       *skEncryptionCKKS
	originalParams          ckks.Parameters // 存储原始CKKS参数

	// 加法全同态加密组件
	additiveHEParams ckks.Parameters // 加法全同态加密参数
	additiveHEPK     *rlwe.PublicKey // 加法全同态加密公钥
	additiveHESK     *rlwe.SecretKey // 加法全同态加密私钥

	// 门限解密组件
	thresholdSKShares map[int]*drlwe.ShamirSecretShare // 门限私钥份额
	thresholdizer     *MyThresholdizer                 // 门限化器

	N int // 总参与方数量
	T int // 门限值
}

// NewFastShamirApproxSS 创建新的FastShamirApproxSS实例
func NewFastShamirApproxSS(N, T int, params ckks.Parameters) (fastShamirApproxSS *FastShamirApproxSS) {
	fastShamirApproxSS = new(FastShamirApproxSS)
	fastShamirApproxSS.N = N
	fastShamirApproxSS.T = T
	fastShamirApproxSS.originalParams = params // 保存原始CKKS参数
	fastShamirApproxSS.VanSS = NewVanillaShamirSS(N, T, params)
	fastShamirApproxSS.params4DoubleEncryption = findLargerParameters(params, N)

	fastShamirApproxSS.f = NewSKencryption(fastShamirApproxSS.params4DoubleEncryption)
	fastShamirApproxSS.params4cmb, _ = rlwe.NewParameters(params.LogN(), params.Q()[:1], nil, 0, params.HammingWeight(), params.Sigma(), params.RingType(), rlwe.NewScale(1), params.DefaultNTTFlag())

	// 初始化加法全同态加密组件
	fastShamirApproxSS.initAdditiveHE()

	// 注意：门限解密组件将在ApproxRecover中初始化，以便准确计算时间

	return
}

// initAdditiveHE 初始化加法全同态加密组件
func (fast *FastShamirApproxSS) initAdditiveHE() {
	// 使用专门为加密原始CKKS私钥设计的参数
	// 而不是用于双重加密的params4DoubleEncryption
	fast.additiveHEParams = generateAdditiveHEParameters(fast.originalParams)

	// 生成加法全同态加密的密钥对
	keyGen := ckks.NewKeyGenerator(fast.additiveHEParams)
	fast.additiveHESK = keyGen.GenSecretKey()
	fast.additiveHEPK = keyGen.GenPublicKey(fast.additiveHESK)
}

// generateThresholdKeys 生成门限解密密钥（使用(T, N) Shamir秘密分享）
func (fast *FastShamirApproxSS) generateThresholdKeys(T, N int) (timeKeyGen time.Duration, sizeCommKey float64) {
	keyGenStart := time.Now()

	fast.thresholdSKShares = make(map[int]*drlwe.ShamirSecretShare)

	// 使用(T, N) Shamir秘密分享
	// 生成N个份额，任意T个可以重构秘密

	points := make([]drlwe.ShamirPublicPoint, N)
	for i := 0; i < N; i++ {
		points[i] = drlwe.ShamirPublicPoint(i + 1)
	}

	// 使用现有的门限化器生成Shamir秘密分享
	fast.thresholdizer = NewMyThresholdizer(fast.additiveHEParams.Parameters)
	shares := fast.thresholdizer.GenShamirSecretShares(T, points, fast.additiveHESK)

	// 存储所有N个份额
	for i := 0; i < N; i++ {
		fast.thresholdSKShares[i] = shares[i]
	}

	timeKeyGen = time.Since(keyGenStart)
	// 通信开销：公钥 + 私钥份额
	sizeCommKey = float64(fast.additiveHEPK.Value[0].MarshalBinarySize64()) +
		float64(fast.thresholdSKShares[0].Poly.Q.MarshalBinarySize64())

	fmt.Printf("Generated (T=%d, N=%d) Shamir secret shares (FastShamirApproxSS)\n", T, N)

	return
}

// additiveHEEncrypt 使用加法全同态加密对密钥进行加密
func (fast *FastShamirApproxSS) additiveHEEncrypt(keyPoly *ring.Poly) *rlwe.Ciphertext {
	encoder := ckks.NewEncoder(fast.additiveHEParams)
	encryptor := ckks.NewEncryptor(fast.additiveHEParams, fast.additiveHEPK)

	// 将密钥多项式转换为复数向量（CKKS需要）
	slots := fast.additiveHEParams.Slots()
	encodeLength := keyPoly.N()
	if encodeLength > slots {
		encodeLength = slots
	}

	keyComplex := make([]complex128, encodeLength)
	for i := 0; i < encodeLength; i++ {
		keyComplex[i] = complex(float64(keyPoly.Coeffs[0][i]), 0)
	}

	// 编码并加密
	plaintext := ckks.NewPlaintext(fast.additiveHEParams, fast.additiveHEParams.MaxLevel())
	encoder.Encode(keyComplex, plaintext, fast.additiveHEParams.LogSlots())

	ciphertext := encryptor.EncryptNew(plaintext)
	return ciphertext
}

// additiveHEAdd 对密文进行同态加法运算
func (fast *FastShamirApproxSS) additiveHEAdd(ct1, ct2 *rlwe.Ciphertext) *rlwe.Ciphertext {
	evaluator := ckks.NewEvaluator(fast.additiveHEParams, rlwe.EvaluationKey{})
	result := ckks.NewCiphertext(fast.additiveHEParams, 1, fast.additiveHEParams.MaxLevel())
	evaluator.Add(ct1, ct2, result)
	return result
}

// additiveHEMulScalar 对密文进行标量乘法运算
func (fast *FastShamirApproxSS) additiveHEMulScalar(ct *rlwe.Ciphertext, scalar uint64) *rlwe.Ciphertext {
	evaluator := ckks.NewEvaluator(fast.additiveHEParams, rlwe.EvaluationKey{})
	encoder := ckks.NewEncoder(fast.additiveHEParams)
	result := ckks.NewCiphertext(fast.additiveHEParams, 1, fast.additiveHEParams.MaxLevel())

	// CKKS需要创建标量明文，然后使用Mul进行乘法
	slots := fast.additiveHEParams.Slots()
	scalarVector := make([]complex128, slots)
	scalarFloat := float64(scalar)
	for i := range scalarVector {
		scalarVector[i] = complex(scalarFloat, 0)
	}
	scalarPt := ckks.NewPlaintext(fast.additiveHEParams, fast.additiveHEParams.MaxLevel())
	encoder.Encode(scalarVector, scalarPt, fast.additiveHEParams.LogSlots())

	// 使用Mul将密文与标量明文相乘
	evaluator.Mul(ct, scalarPt, result)
	return result
}

// thresholdDecryptPartial 计算部分解密份额（先乘以拉格朗日系数，再添加smudging noise）
func (fast *FastShamirApproxSS) thresholdDecryptPartial(participantID int, ciphertext *rlwe.Ciphertext, lagrangeCoeff uint64, bound4smudgingNoise int) (*rlwe.Ciphertext, time.Duration) {
	decStart := time.Now()

	// 使用参与方的私钥份额进行部分解密
	share := fast.thresholdSKShares[participantID]

	// 创建部分解密结果
	partialDecryption := ckks.NewCiphertext(fast.additiveHEParams, 0, fast.additiveHEParams.MaxLevel())

	// 计算 c0 + share.sk * c1
	ringQ := fast.additiveHEParams.RingQ()
	ringQ.MulCoeffsMontgomeryLvl(ciphertext.Level(), ciphertext.Value[1], share.Poly.Q, partialDecryption.Value[0])
	ringQ.Add(ciphertext.Value[0], partialDecryption.Value[0], partialDecryption.Value[0])

	// 先乘以拉格朗日系数
	lambdaScalar := ringQ.NewRNSScalarFromUInt64(lagrangeCoeff)
	temp := ringQ.NewPoly()
	ringQ.MulRNSScalarMontgomery(partialDecryption.Value[0], lambdaScalar, temp)
	partialDecryption.Value[0] = temp

	// 添加smudging noise以保护隐私
	smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSampler(fast.originalParams, ringQ, bound4smudgingNoise)
	noise := ringQ.NewPoly()
	smudgingNoiseSampler.ReadLvl(fast.additiveHEParams.MaxLevel(), noise)

	if actualBound != bound4smudgingNoise {
		fmt.Printf("FastShamirApproxSS thresholdDecryptPartial: Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
	}

	// 将noise添加到部分解密结果中
	ringQ.Add(partialDecryption.Value[0], noise, partialDecryption.Value[0])

	decTime := time.Since(decStart)
	return partialDecryption, decTime
}

// thresholdDecryptCombine 组合部分解密份额得到最终结果（Shamir秘密分享版本）
func (fast *FastShamirApproxSS) thresholdDecryptCombine(participants []int, partialDecryptions map[int]*rlwe.Ciphertext, bound4smudgingNoise int) *ring.Poly {
	// 对于Shamir秘密分享的T-out-of-N门限解密，直接将所有部分解密结果相加
	// 注意：每个参与方已经在自己的部分解密中乘以了拉格朗日系数并添加了smudging noise

	ringQ := fast.additiveHEParams.RingQ()
	result := ringQ.NewPoly()
	result.Zero()

	// 直接将所有带噪声的部分解密结果相加
	for _, p := range participants {
		// 部分解密结果已经乘以拉格朗日系数并包含了噪声，直接相加
		ringQ.Add(result, partialDecryptions[p].Value[0], result)
	}

	return result
}

// ApproxRecover 基于加法全同态加密和门限解密的近似恢复算法
func (fast *FastShamirApproxSS) ApproxRecover(bound4smudgingNoise int) (isSuc bool, timeCompOnce time.Duration, timeCompNotOnce time.Duration, sizeCommOnce float64, sizeCommNotOnce float64) {

	paramsMesSpace := fast.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
	var timeNotOnce1, timeNotOnce2 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64
	var sizeCommNotOnce1 float64

	//Round 1 Starts - 使用加法全同态加密替代秘密分享

	//Samples the polynomial a in our paper
	CRS, _ := utils.NewPRNG()
	a := TagGen(fast.params4DoubleEncryption, CRS)
	//Select the participants in round 1
	parR1 := SelectRandomParticipants(fast.VanSS.N, fast.VanSS.T)
	//The output of parties and aggregator in Round 1
	CTni_all := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)
	CTsi_all := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)
	lagrangeCoeffs := make(map[int]uint64, fast.VanSS.T)

	// 存储密钥的密文
	ek1Ciphertexts := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)
	ek2Ciphertexts := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)

	//Parties run PartyR1 - 使用加法全同态加密
	for _, par := range parR1 {
		//Read the share to be encrypted from file
		filename := fmt.Sprintf("%s%d", "mesToParty", par)
		skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
		skFile, _ := os.Open(skFilename)
		defer skFile.Close()
		skData, _ := io.ReadAll(skFile)
		share := fast.VanSS.thdizer.AllocateThresholdSecretShare()
		share.UnmarshalBinary(skData)

		// 生成密钥但不进行秘密分享，而是进行加法全同态加密
		timeStage1Start := time.Now()
		ek1, _ := fast.f.GenerateEKOnly()
		ek2, _ := fast.f.GenerateEKOnly()

		// 对密钥进行加法全同态加密
		ek1Ciphertexts[par] = fast.additiveHEEncrypt(ek1.Value.Q)
		ek2Ciphertexts[par] = fast.additiveHEEncrypt(ek2.Value.Q)

		timeStage1 = time.Since(timeStage1Start)
		sizeCommStage1 = float64(ek1Ciphertexts[par].Value[0].MarshalBinarySize64() * 2) // ek1和ek2的密文大小

		timeNotOnce1Start := time.Now()
		// 使用安全的smudging noise采样器
		smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSampler(fast.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)

		ni := ringQ.NewPoly()
		smudgingNoiseSampler.ReadLvl(levelQ, ni)

		// 输出安全信息
		if actualBound != bound4smudgingNoise {
			fmt.Printf("FastShamirApproxSS: Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
		}

		CTsi_all[par] = fast.f.EncPolyCoeff(a, ek1, share.Poly.Q)
		CTni_all[par] = fast.f.EncPolyCoeff(a, ek2, ni)
		timeNotOnce1 = time.Since(timeNotOnce1Start)
		sizeCommNotOnce1 = float64(CTsi_all[par].Value[0].MarshalBinarySize64()) * 2
	}

	//AggregatorR1 - 计算拉格朗日系数
	timeStage2Start := time.Now()
	survivialPublicPoint := []drlwe.ShamirPublicPoint{}
	//Collect survivialClient
	for _, pi := range parR1 {
		survivialPublicPoint = append(survivialPublicPoint, drlwe.ShamirPublicPoint(pi+1))
	}
	cmbR1 := NewMyCombiner(&fast.params4cmb, survivialPublicPoint, fast.VanSS.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(fast.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}
	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(8 * len(lagrangeCoeffs))
	//Round 1 Ends!

	//Round 2 Starts - 使用门限解密替代原始的DK生成

	// 生成(T, N)门限解密密钥
	timeKeyGen, sizeCommKey := fast.generateThresholdKeys(fast.T, fast.N)

	// 不再选择新的参与者，直接使用parR1集合进行T-out-of-N门限解密
	// 所有parR1中的T个参与者都必须参与门限解密过程

	timeStage3Start := time.Now()
	// 计算 ∑_{i∈T} λ_i * ek_{i,1} + ek_{i,2} 的密文
	combinedCiphertext := fast.additiveHEEncrypt(fast.additiveHEParams.RingQ().NewPoly()) // 初始化为0的密文
	combinedCiphertext.Value[0].Zero()
	combinedCiphertext.Value[1].Zero()

	for _, cIdx := range parR1 {
		// λ_i * ek_{i,1}
		scaledEk1 := fast.additiveHEMulScalar(ek1Ciphertexts[cIdx], lagrangeCoeffs[cIdx])
		combinedCiphertext = fast.additiveHEAdd(combinedCiphertext, scaledEk1)

		// + ek_{i,2}
		combinedCiphertext = fast.additiveHEAdd(combinedCiphertext, ek2Ciphertexts[cIdx])
	}
	timeStage3 = time.Since(timeStage3Start)
	sizeCommStage3 = float64(combinedCiphertext.Value[0].MarshalBinarySize64())

	// T-out-of-N门限解密阶段，所有parR1中的参与者都参与
	timeStage4Start := time.Now()
	partialDecryptions := make(map[int]*rlwe.Ciphertext)

	// 建立parR1参与者到门限私钥份额的映射
	// parR1中第i个参与者使用门限私钥份额parR1[i]
	for i, parID := range parR1 {
		// 获取对应的拉格朗日系数
		lagrangeCoeff := lagrangeCoeffs[parID]
		partialDecryptions[i], _ = fast.thresholdDecryptPartial(parID, combinedCiphertext, lagrangeCoeff, bound4smudgingNoise)
	}

	// 组合部分解密结果得到 ∑_{i∈T} λ_i * ek_{i,1} + ek_{i,2}
	// 使用连续的ID [0, 1, 2, ..., T-1] 进行组合
	thresholdParticipants := make([]int, fast.T)
	for i := 0; i < fast.T; i++ {
		thresholdParticipants[i] = i
	}
	combinedKey := fast.thresholdDecryptCombine(thresholdParticipants, partialDecryptions, bound4smudgingNoise)
	timeStage4 = time.Since(timeStage4Start)

	// 使用组合后的密钥进行最终解密
	timeNotOnce2Start := time.Now()
	approxMessage := ringQ.NewPoly()
	// 使用专门的加法全同态解密函数
	fast.f.FEDecFinalWithAdditiveHEKey(combinedKey, fast.additiveHEParams, CTni_all, CTsi_all, lagrangeCoeffs, parR1, approxMessage)
	timeNotOnce2 = time.Since(timeNotOnce2Start)

	//Round 2 Ends!

	timeCompOnce = timeStage1 + timeStage2 + timeKeyGen + timeStage3 + timeStage4
	timeCompNotOnce = timeNotOnce1 + timeNotOnce2
	sizeCommOnce = sizeCommStage1 + sizeCommStage2 + sizeCommKey + sizeCommStage3
	sizeCommNotOnce = sizeCommNotOnce1
	isSuc = TestResult(approxMessage, fast.VanSS.secret.Value.Q, uint64(fast.VanSS.T*bound4smudgingNoise), fast.VanSS.thdizer.params.Q()[0:1])
	return
}

// ApproxRecover4TestTime 用于测试时间的版本，优化了重复操作
func (fast *FastShamirApproxSS) ApproxRecover4TestTime(bound4smudgingNoise int) (timeCompOnce time.Duration, timeCompNotOnce time.Duration, sizeCommOnce float64, sizeCommNotOnce float64) {

	paramsMesSpace := fast.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
	var timeNotOnce1, timeNotOnce2 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64
	var sizeCommNotOnce1 float64

	//Round 1 Starts

	//Samples the polynomial a in our paper
	CRS, _ := utils.NewPRNG()
	a := TagGen(fast.params4DoubleEncryption, CRS)
	//Select the participants in round 1
	parR1 := SelectRandomParticipants(fast.VanSS.N, fast.VanSS.T)
	//The output of parties and aggregator in Round 1
	CTni_all := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)
	CTsi_all := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)
	lagrangeCoeffs := make(map[int]uint64, fast.VanSS.T)

	// 存储加密后的密钥密文
	ek1Ciphertexts := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)
	ek2Ciphertexts := make(map[int]*rlwe.Ciphertext, fast.VanSS.T)

	//Parties run PartyR1 - 测试时间版本，只计算第一个参与方的时间
	for i, par := range parR1 {
		if i == 0 {
			//Read the share to be encrypted from file
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			defer skFile.Close()
			skData, _ := io.ReadAll(skFile)
			share := fast.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			timeStage1Start := time.Now()
			// 生成ek1和ek2
			keyGen := ckks.NewKeyGenerator(fast.params4DoubleEncryption)
			ek1Orig := keyGen.GenSecretKey()
			ek2Orig := keyGen.GenSecretKey()

			// 转换为a·ek+e'形式
			smudgingNoiseSampler, _ := CreateSecureSmudgingNoiseSampler(fast.originalParams, fast.params4DoubleEncryption.RingQ(), bound4smudgingNoise)

			// 对ek1: a·ek1 + e1'
			e1Prime := fast.params4DoubleEncryption.RingQ().NewPoly()
			smudgingNoiseSampler.ReadLvl(fast.params4DoubleEncryption.MaxLevel(), e1Prime)
			convertedEk1 := fast.params4DoubleEncryption.RingQ().NewPoly()
			fast.params4DoubleEncryption.RingQ().MulCoeffsMontgomeryLvl(fast.params4DoubleEncryption.MaxLevel(), a, ek1Orig.Value.Q, convertedEk1)
			fast.params4DoubleEncryption.RingQ().AddLvl(fast.params4DoubleEncryption.MaxLevel(), convertedEk1, e1Prime, convertedEk1)

			// 对ek2: a·ek2 + e2'
			e2Prime := fast.params4DoubleEncryption.RingQ().NewPoly()
			smudgingNoiseSampler.ReadLvl(fast.params4DoubleEncryption.MaxLevel(), e2Prime)
			convertedEk2 := fast.params4DoubleEncryption.RingQ().NewPoly()
			fast.params4DoubleEncryption.RingQ().MulCoeffsMontgomeryLvl(fast.params4DoubleEncryption.MaxLevel(), a, ek2Orig.Value.Q, convertedEk2)
			fast.params4DoubleEncryption.RingQ().AddLvl(fast.params4DoubleEncryption.MaxLevel(), convertedEk2, e2Prime, convertedEk2)

			// 对转换后的份额进行加密（使用加法全同态加密）
			ek1Ciphertexts[par] = fast.additiveHEEncrypt(convertedEk1)
			ek2Ciphertexts[par] = fast.additiveHEEncrypt(convertedEk2)

			timeStage1 = time.Since(timeStage1Start)
			// 计算密文的通信大小
			sizeCommStage1 = float64(ek1Ciphertexts[par].Value[0].MarshalBinarySize64()*2) * 2 // ek1和ek2的密文

			timeNotOnce1Start := time.Now()
			// 使用安全的smudging noise采样器
			smudgingNoiseSampler2, actualBound := CreateSecureSmudgingNoiseSampler(fast.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)

			ni := ringQ.NewPoly()
			smudgingNoiseSampler2.ReadLvl(levelQ, ni)

			// 输出安全信息
			if actualBound != bound4smudgingNoise {
				fmt.Printf("FastShamirApproxSS (TestTime): Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
			}

			CTsi_all[par] = fast.f.EncPolyCoeff(a, ek1Orig, share.Poly.Q)
			CTni_all[par] = fast.f.EncPolyCoeff(a, ek2Orig, ni)
			timeNotOnce1 = time.Since(timeNotOnce1Start)
			sizeCommNotOnce1 = float64(CTsi_all[par].Value[0].MarshalBinarySize64()) * 2
		} else {
			// 复用第一个参与方的结果以节省测试时间
			ek1Ciphertexts[par] = ek1Ciphertexts[parR1[0]]
			ek2Ciphertexts[par] = ek2Ciphertexts[parR1[0]]
			CTsi_all[par] = CTsi_all[parR1[0]]
			CTni_all[par] = CTni_all[parR1[0]]
		}
	}

	//AggregatorR1
	timeStage2Start := time.Now()
	survivialPublicPoint := []drlwe.ShamirPublicPoint{}
	for _, pi := range parR1 {
		survivialPublicPoint = append(survivialPublicPoint, drlwe.ShamirPublicPoint(pi+1))
	}
	cmbR1 := NewMyCombiner(&fast.params4cmb, survivialPublicPoint, fast.VanSS.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(fast.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}
	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(8 * len(lagrangeCoeffs))
	//Round 1 Ends!

	//Round 2 Starts!

	// 生成(T, N)门限解密密钥
	timeKeyGen, sizeCommKey := fast.generateThresholdKeys(fast.T, fast.N)

	// 对密文进行同态运算：计算 ∑_{i∈T} λ_i * ek_{i,1} + ek_{i,2} 的密文
	// 只对第一个参与方记录时间
	timeStage3Start := time.Now()
	combinedCiphertext := fast.additiveHEEncrypt(fast.additiveHEParams.RingQ().NewPoly()) // 初始化为0的密文
	combinedCiphertext.Value[0].Zero()
	combinedCiphertext.Value[1].Zero()

	for i, cIdx := range parR1 {
		if i == 0 {
			// λ_i * ek_{i,1}
			scaledEk1 := fast.additiveHEMulScalar(ek1Ciphertexts[cIdx], lagrangeCoeffs[cIdx])
			combinedCiphertext = fast.additiveHEAdd(combinedCiphertext, scaledEk1)

			// + ek_{i,2}
			combinedCiphertext = fast.additiveHEAdd(combinedCiphertext, ek2Ciphertexts[cIdx])
			timeStage3 = time.Since(timeStage3Start)
		} else {
			// 其他参与方复用计算
			scaledEk1 := fast.additiveHEMulScalar(ek1Ciphertexts[cIdx], lagrangeCoeffs[cIdx])
			combinedCiphertext = fast.additiveHEAdd(combinedCiphertext, scaledEk1)
			combinedCiphertext = fast.additiveHEAdd(combinedCiphertext, ek2Ciphertexts[cIdx])
		}
	}
	sizeCommStage3 = float64(combinedCiphertext.Value[0].MarshalBinarySize64())

	// (T, N)门限解密阶段
	timeStage4Start := time.Now()
	partialDecryptions := make(map[int]*rlwe.Ciphertext)

	// 只对第一个参与方执行真正的部分解密，其他参与方复用结果
	// 建立parR1参与者到门限私钥份额的映射
	for i, parID := range parR1 {
		if i == 0 {
			// 获取对应的拉格朗日系数
			lagrangeCoeff := lagrangeCoeffs[parID]
			partialDecryptions[i], _ = fast.thresholdDecryptPartial(parID, combinedCiphertext, lagrangeCoeff, bound4smudgingNoise)
		} else {
			partialDecryptions[i] = partialDecryptions[0]
		}
	}

	// 组合部分解密结果得到总私钥
	thresholdParticipants := make([]int, fast.T)
	for i := 0; i < fast.T; i++ {
		thresholdParticipants[i] = i
	}
	combinedKey := fast.thresholdDecryptCombine(thresholdParticipants, partialDecryptions, bound4smudgingNoise)
	timeStage4 = time.Since(timeStage4Start)

	timeNotOnce2Start := time.Now()
	approxMessage := ringQ.NewPoly()
	// 使用组合后的密钥进行最终解密
	// combinedKey 需要从 additiveHEParams 的环映射到 params4DoubleEncryption 的环
	fast.f.FEDecFinalWithAdditiveHEKey(combinedKey, fast.additiveHEParams, CTni_all, CTsi_all, lagrangeCoeffs, parR1, approxMessage)
	timeNotOnce2 = time.Since(timeNotOnce2Start)

	//Round 2 Ends!
	//timeKeyGen可以预处理阶段进行
	timeCompOnce = timeStage1 + timeStage2 + timeStage3 + timeStage4
	timeCompNotOnce = timeNotOnce1 + timeNotOnce2
	sizeCommOnce = sizeCommStage1 + sizeCommStage2 + sizeCommKey + sizeCommStage3
	sizeCommNotOnce = sizeCommNotOnce1

	fmt.Println("timeStageFast：", timeStage1, timeStage2, timeKeyGen, timeStage3, timeStage4)
	fmt.Println("timeNotOnceFast： ", timeNotOnce1, timeNotOnce2)
	return
}
