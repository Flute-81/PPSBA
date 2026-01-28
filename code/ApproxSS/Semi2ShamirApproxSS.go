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

// Semi2ShamirApproxSS 基于(T,T)门限全同态加密的半门限Shamir近似秘密共享方案
type Semi2ShamirApproxSS struct {
	VanSS *VanillaShamirSS

	params4cmb        rlwe.Parameters
	thresholdizer     *MyThresholdizer           // 门限化器
	thresholdSKShares []*drlwe.ShamirSecretShare // (T,T)门限私钥份额
	publicKey         *rlwe.PublicKey            // 门限公钥
	N                 int                        // 总参与方数量
	T                 int                        // 门限值
}

// NewSemi2ShamirApproxSS 创建新的Semi2ShamirApproxSS实例
func NewSemi2ShamirApproxSS(N, T int, params ckks.Parameters) (semi2ShamirApproxSS *Semi2ShamirApproxSS) {
	semi2ShamirApproxSS = new(Semi2ShamirApproxSS)
	semi2ShamirApproxSS.N = N
	semi2ShamirApproxSS.T = T
	semi2ShamirApproxSS.VanSS = NewVanillaShamirSS(N, T, params)

	// 设置组合器参数，使用第一个模数创建新的参数
	semi2ShamirApproxSS.params4cmb, _ = rlwe.NewParameters(params.LogN(), params.Q()[:1], nil, 0, params.HammingWeight(), params.Sigma(), params.RingType(), rlwe.NewScale(1), params.DefaultNTTFlag())

	// 初始化门限化器
	semi2ShamirApproxSS.thresholdizer = NewMyThresholdizer(semi2ShamirApproxSS.params4cmb)

	return
}

// generateThresholdKeys 生成(T,T)门限密钥对（使用加性秘密分享）
func (semi2 *Semi2ShamirApproxSS) generateThresholdKeys() (timeKeyGen time.Duration, sizeCommKey float64) {
	keyGenStart := time.Now()

	// 1. 生成主私钥
	kgen := rlwe.NewKeyGenerator(semi2.params4cmb)
	sk := kgen.GenSecretKey()

	// 2. 生成公钥
	semi2.publicKey = kgen.GenPublicKey(sk)

	// 3. 使用加性秘密分享私钥：sk = s1 + s2 + ... + sT (mod q)
	ringQ := semi2.params4cmb.RingQ()
	semi2.thresholdSKShares = make([]*drlwe.ShamirSecretShare, semi2.T)

	// 生成前T-1个随机份额（与FastShamirApproxSS完全相同的策略）
	prng, _ := utils.NewPRNG()
	gaussianSampler := ring.NewGaussianSampler(prng, ringQ, semi2.params4cmb.Sigma(), int(6*semi2.params4cmb.Sigma()))

	shares := make([]*ring.Poly, semi2.T)
	sum := ringQ.NewPoly()
	sum.Zero()

	for i := 0; i < semi2.T-1; i++ {
		shares[i] = ringQ.NewPoly()
		// 使用与FastShamirApproxSS完全相同的采样方式
		gaussianSampler.ReadLvl(semi2.params4cmb.MaxLevel(), shares[i])
		// 累加到sum中
		ringQ.Add(sum, shares[i], sum)

		// 转换为ShamirSecretShare格式以保持接口一致
		semi2.thresholdSKShares[i] = semi2.thresholdizer.AllocateThresholdSecretShare()
		semi2.thresholdSKShares[i].Q = shares[i].CopyNew()
	}

	// 第T个份额 = sk - (s1 + s2 + ... + s_{T-1}) (mod q)
	semi2.thresholdSKShares[semi2.T-1] = semi2.thresholdizer.AllocateThresholdSecretShare()
	shares[semi2.T-1] = ringQ.NewPoly()
	ringQ.Sub(sk.Value.Q, sum, shares[semi2.T-1])
	semi2.thresholdSKShares[semi2.T-1].Q.Copy(shares[semi2.T-1])

	timeKeyGen = time.Since(keyGenStart)
	// 通信开销：公钥 + T个私钥份额
	sizeCommKey = float64(semi2.publicKey.Value[0].MarshalBinarySize64()) +
		float64(semi2.thresholdSKShares[0].Q.MarshalBinarySize64())

	fmt.Printf("Generated additive secret shares for %d participants\n", semi2.T)

	return
}

// thresholdEncrypt 使用门限公钥进行加密
func (semi2 *Semi2ShamirApproxSS) thresholdEncrypt(plaintext *ring.Poly) (*rlwe.Ciphertext, time.Duration) {
	encStart := time.Now()

	encryptor := rlwe.NewEncryptor(semi2.params4cmb, semi2.publicKey)
	pt := rlwe.NewPlaintext(semi2.params4cmb, semi2.params4cmb.MaxLevel())
	pt.Value.Copy(plaintext)

	ciphertext := encryptor.EncryptNew(pt)

	encTime := time.Since(encStart)
	return ciphertext, encTime
}

// homomorphicAdd 同态加法
func (semi2 *Semi2ShamirApproxSS) homomorphicAdd(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, time.Duration) {
	addStart := time.Now()

	result := ct1.CopyNew()
	ringQ := semi2.params4cmb.RingQ()
	ringQ.Add(result.Value[0], ct2.Value[0], result.Value[0])
	ringQ.Add(result.Value[1], ct2.Value[1], result.Value[1])

	addTime := time.Since(addStart)
	return result, addTime
}

// homomorphicMulScalar 同态标量乘法
func (semi2 *Semi2ShamirApproxSS) homomorphicMulScalar(ct *rlwe.Ciphertext, scalar uint64) (*rlwe.Ciphertext, time.Duration) {
	mulStart := time.Now()

	result := ct.CopyNew()
	ringQ := semi2.params4cmb.RingQ()
	ringQ.MulScalar(result.Value[0], scalar, result.Value[0])
	ringQ.MulScalar(result.Value[1], scalar, result.Value[1])

	mulTime := time.Since(mulStart)
	return result, mulTime
}

// thresholdDecryptPartial 门限部分解密
func (semi2 *Semi2ShamirApproxSS) thresholdDecryptPartial(participantID int, ciphertext *rlwe.Ciphertext, bound4smudgingNoise int) (*rlwe.Ciphertext, time.Duration) {
	decStart := time.Now()

	// 使用参与方的私钥份额进行部分解密
	share := semi2.thresholdSKShares[participantID]

	// 创建部分解密结果
	partialDecryption := rlwe.NewCiphertext(semi2.params4cmb, 0, semi2.params4cmb.MaxLevel())

	// 计算 c0 + share.sk * c1
	ringQ := semi2.params4cmb.RingQ()
	ringQ.MulCoeffsMontgomeryLvl(ciphertext.Level(), ciphertext.Value[1], share.Q, partialDecryption.Value[0])
	ringQ.Add(ciphertext.Value[0], partialDecryption.Value[0], partialDecryption.Value[0])

	// 添加smudging noise以保护隐私
	smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSamplerFromRLWE(&semi2.params4cmb, ringQ, bound4smudgingNoise)
	noise := ringQ.NewPoly()
	smudgingNoiseSampler.ReadLvl(semi2.params4cmb.MaxLevel(), noise)

	if actualBound != bound4smudgingNoise {
		fmt.Printf("Semi2ShamirApproxSS thresholdDecryptPartial: Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
	}

	// 将noise添加到部分解密结果中
	ringQ.Add(partialDecryption.Value[0], noise, partialDecryption.Value[0])

	decTime := time.Since(decStart)
	return partialDecryption, decTime
}

// thresholdDecryptCombine 组合门限部分解密结果（加性秘密分享版本）
func (semi2 *Semi2ShamirApproxSS) thresholdDecryptCombine(partialDecryptions []*rlwe.Ciphertext) (*ring.Poly, time.Duration) {
	combineStart := time.Now()

	// 对于加性秘密分享的(T,T)门限解密，直接将所有部分解密结果相加
	// 因为 sk = s1 + s2 + ... + sT，所以解密结果也是相加关系
	ringQ := semi2.params4cmb.RingQ()
	result := ringQ.NewPoly()
	result.Zero()

	for _, partialDec := range partialDecryptions {
		ringQ.Add(result, partialDec.Value[0], result)
	}

	combineTime := time.Since(combineStart)
	return result, combineTime
}

// ApproxRecover 执行三轮门限全同态加密近似恢复算法
func (semi2 *Semi2ShamirApproxSS) ApproxRecover(bound4smudgingNoise int) (timeComp time.Duration, sizeComm float64) {

	paramsMesSpace := semi2.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64

	// Round 1: 参与方发送标记信息，聚合者计算拉格朗日系数
	fmt.Println("=== Round 1: 标记信息收集与拉格朗日系数计算 ===")

	// 选择T个参与方
	parR1 := SelectRandomParticipants(semi2.N, semi2.T)
	fmt.Printf("Selected participants: %v\n", parR1)

	timeStage1Start := time.Now()

	// 聚合者收集所有标记信息，计算拉格朗日系数
	survivialPublicPoint := make([]drlwe.ShamirPublicPoint, len(parR1))
	for i, pi := range parR1 {
		survivialPublicPoint[i] = drlwe.ShamirPublicPoint(pi + 1)
	}

	// 计算拉格朗日系数
	cmbR1 := NewMyCombiner(&semi2.params4cmb, survivialPublicPoint, semi2.T)
	lagrangeCoeffs := make(map[int]uint64, semi2.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(semi2.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}

	timeStage1 = time.Since(timeStage1Start)
	sizeCommStage1 = float64(4) + float64(8*len(lagrangeCoeffs)) // 参与方ID + 拉格朗日系数

	fmt.Printf("Round 1 completed: Lagrange coefficients computed for %d participants\n", len(lagrangeCoeffs))

	// Round 2: (T,T)门限全同态加密
	fmt.Println("=== Round 2: (T,T)门限全同态加密 ===")

	timeStage2Start := time.Now()

	// 生成(T,T)门限密钥对
	timeKeyGen, sizeCommKey := semi2.generateThresholdKeys()

	// 存储加密后的数据
	encryptedBi := make(map[int]*rlwe.Ciphertext, semi2.T)
	encryptedNi := make(map[int]*rlwe.Ciphertext, semi2.T)

	// 记录单个参与方的加密时间
	var singleEncryptTime time.Duration

	for _, par := range parR1 {
		// 读取参与方的秘密份额 bᵢ
		filename := fmt.Sprintf("%s%d", "mesToParty", par)
		skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
		skFile, _ := os.Open(skFilename)
		defer skFile.Close()
		skData, _ := io.ReadAll(skFile)
		share := semi2.VanSS.thdizer.AllocateThresholdSecretShare()
		share.UnmarshalBinary(skData)

		// 生成噪声 nᵢ
		smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSamplerFromRLWE(semi2.VanSS.thdizer.params, paramsMesSpace.RingQ(), bound4smudgingNoise)
		ni := ringQ.NewPoly()
		smudgingNoiseSampler.ReadLvl(levelQ, ni)

		if actualBound != bound4smudgingNoise {
			fmt.Printf("Semi2ShamirApproxSS: Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
		}

		// 对 bᵢ 和 nᵢ 进行门限全同态加密
		var encTimeB, encTimeN time.Duration
		encryptedBi[par], encTimeB = semi2.thresholdEncrypt(share.Poly.Q)
		encryptedNi[par], encTimeN = semi2.thresholdEncrypt(ni)

		// 记录加密时间（只记录第一个参与方的时间）
		if par == parR1[0] {
			singleEncryptTime = encTimeB + encTimeN
		}

		fmt.Printf("Participant %d: encrypted bᵢ and nᵢ\n", par)
	}

	// 聚合者同态计算 ∑_{i∈T} λᵢ^(T)·bᵢ + nᵢ
	var combinedCiphertext *rlwe.Ciphertext
	var singleHomomorphicTime time.Duration

	for idx, par := range parR1 {
		lagrangeCoeff := lagrangeCoeffs[par]

		// 计算 λᵢ·bᵢ
		scaledBi, mulTime := semi2.homomorphicMulScalar(encryptedBi[par], lagrangeCoeff)

		// 计算 λᵢ·bᵢ + nᵢ
		partialResult, addTime1 := semi2.homomorphicAdd(scaledBi, encryptedNi[par])

		if idx == 0 {
			combinedCiphertext = partialResult
			singleHomomorphicTime = mulTime + addTime1
		} else {
			var addTime2 time.Duration
			combinedCiphertext, addTime2 = semi2.homomorphicAdd(combinedCiphertext, partialResult)
			if idx == 1 {
				singleHomomorphicTime += addTime2
			}
		}

		fmt.Printf("Participant %d: added λ%d·b%d + n%d to combined result\n", par, par, par, par)
	}

	timeStage2 = time.Since(timeStage2Start)
	// 通信开销：密钥生成 + 加密数据
	sizeCommStage2 = sizeCommKey + float64(encryptedBi[parR1[0]].Value[0].MarshalBinarySize64())*2

	fmt.Printf("Round 2 completed: Combined ciphertext computed\n")

	// Round 3: 门限解密
	fmt.Println("=== Round 3: 门限解密 ===")

	timeStage3Start := time.Now()

	// T个参与方进行门限部分解密
	partialDecryptions := make([]*rlwe.Ciphertext, semi2.T)
	var singlePartialDecTime time.Duration

	for i := 0; i < semi2.T; i++ {
		partialDecryptions[i], singlePartialDecTime = semi2.thresholdDecryptPartial(i, combinedCiphertext, bound4smudgingNoise)
		fmt.Printf("Participant %d: computed partial decryption\n", i)
	}

	// 聚合者组合部分解密结果
	approxMessage, combineTime := semi2.thresholdDecryptCombine(partialDecryptions)

	timeStage3 = time.Since(timeStage3Start)
	// 通信开销：T个部分解密结果
	sizeCommStage3 = float64(partialDecryptions[0].Value[0].MarshalBinarySize64())

	fmt.Printf("Round 3 completed: Final result decrypted\n")

	// 计算总时间和通信量
	timeComp = timeStage1 + timeStage2 + timeStage3
	sizeComm = sizeCommStage1 + sizeCommStage2 + sizeCommStage3

	fmt.Printf("=== 性能统计 ===\n")
	fmt.Printf("Time breakdown:\n")
	fmt.Printf("  Stage1 (标记&拉格朗日): %v\n", timeStage1)
	fmt.Printf("  Stage2 (密钥生成=%v, 加密=%v, 同态计算=%v): %v\n", timeKeyGen, singleEncryptTime, singleHomomorphicTime, timeStage2)
	fmt.Printf("  Stage3 (部分解密=%v, 组合=%v): %v\n", singlePartialDecTime, combineTime, timeStage3)
	fmt.Printf("  Total: %v\n", timeComp)
	fmt.Printf("Communication breakdown:\n")
	fmt.Printf("  Stage1: %.2f bytes\n", sizeCommStage1)
	fmt.Printf("  Stage2: %.2f bytes (keys=%.2f, encryptions=%.2f)\n", sizeCommStage2, sizeCommKey, sizeCommStage2-sizeCommKey)
	fmt.Printf("  Stage3: %.2f bytes\n", sizeCommStage3)
	fmt.Printf("  Total: %.2f bytes\n", sizeComm)

	// 测试恢复结果
	_ = TestResult(approxMessage, semi2.VanSS.secret.Value.Q, uint64(semi2.T*bound4smudgingNoise), semi2.VanSS.thdizer.params.Q()[0:1])

	return
}

// ApproxRecover4TestTime 用于测试时间的版本（只处理第一个参与方的实际计算）
func (semi2 *Semi2ShamirApproxSS) ApproxRecover4TestTime(bound4smudgingNoise int) (timeComp time.Duration, sizeComm float64) {

	paramsMesSpace := semi2.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64

	// Round 1: 参与方发送标记信息，聚合者计算拉格朗日系数
	parR1 := SelectRandomParticipants(semi2.N, semi2.T)

	timeStage1Start := time.Now()

	survivialPublicPoint := make([]drlwe.ShamirPublicPoint, len(parR1))
	for i, pi := range parR1 {
		survivialPublicPoint[i] = drlwe.ShamirPublicPoint(pi + 1)
	}

	cmbR1 := NewMyCombiner(&semi2.params4cmb, survivialPublicPoint, semi2.T)
	lagrangeCoeffs := make(map[int]uint64, semi2.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(semi2.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}

	timeStage1 = time.Since(timeStage1Start)
	sizeCommStage1 = float64(4) + float64(8*len(lagrangeCoeffs))

	// Round 2: (T,T)门限全同态加密（测试版本）
	timeStage2Start := time.Now()

	timeKeyGen, sizeCommKey := semi2.generateThresholdKeys()

	encryptedBi := make(map[int]*rlwe.Ciphertext, semi2.T)
	encryptedNi := make(map[int]*rlwe.Ciphertext, semi2.T)

	// 只对第一个参与方进行实际加密
	for i, par := range parR1 {
		if i == 0 {
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			defer skFile.Close()
			skData, _ := io.ReadAll(skFile)
			share := semi2.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSamplerFromRLWE(semi2.VanSS.thdizer.params, paramsMesSpace.RingQ(), bound4smudgingNoise)
			ni := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, ni)

			if actualBound != bound4smudgingNoise {
				fmt.Printf("Semi2ShamirApproxSS (TestTime): Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
			}

			encryptedBi[par], _ = semi2.thresholdEncrypt(share.Poly.Q)
			encryptedNi[par], _ = semi2.thresholdEncrypt(ni)
		} else {
			// 其他参与方复用第一个参与方的结果
			encryptedBi[par] = encryptedBi[parR1[0]]
			encryptedNi[par] = encryptedNi[parR1[0]]
		}
	}

	// 同态计算（只计算第一个参与方的实际时间）
	lagrangeCoeff := lagrangeCoeffs[parR1[0]]
	scaledBi, _ := semi2.homomorphicMulScalar(encryptedBi[parR1[0]], lagrangeCoeff)
	combinedCiphertext, _ := semi2.homomorphicAdd(scaledBi, encryptedNi[parR1[0]])

	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = sizeCommKey + float64(encryptedBi[parR1[0]].Value[0].MarshalBinarySize64())*2

	// Round 3: 门限解密（测试版本）
	timeStage3Start := time.Now()

	// 只对第一个参与方进行实际部分解密，其他参与方复用结果
	partialDecryptions := make([]*rlwe.Ciphertext, semi2.T)
	partialDecryptions[0], _ = semi2.thresholdDecryptPartial(0, combinedCiphertext, bound4smudgingNoise)

	// 其他参与方复用结果
	for i := 1; i < semi2.T; i++ {
		partialDecryptions[i] = partialDecryptions[0]
	}

	// 组合部分解密结果
	_, _ = semi2.thresholdDecryptCombine(partialDecryptions)

	timeStage3 = time.Since(timeStage3Start)
	sizeCommStage3 = float64(partialDecryptions[0].Value[0].MarshalBinarySize64())

	// 计算总时间和通信量
	timeComp = timeStage1 + timeStage2 + timeStage3
	sizeComm = sizeCommStage1 + sizeCommStage2 + sizeCommStage3

	fmt.Printf("TestTime - Time breakdown: Stage1=%v, Stage2=%v(KeyGen=%v), Stage3=%v, Total=%v\n",
		timeStage1, timeStage2, timeKeyGen, timeStage3, timeComp)
	fmt.Printf("TestTime - Size breakdown: Stage1=%v, Stage2=%v, Stage3=%v, Total=%v\n",
		sizeCommStage1, sizeCommStage2, sizeCommStage3, sizeComm)

	return
}
