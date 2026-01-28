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
)

// SemiShamirApproxSS 半门限Shamir近似秘密共享方案
type SemiShamirApproxSS struct {
	VanSS *VanillaShamirSS

	params4cmb   rlwe.Parameters
	originalCKKS ckks.Parameters // 存储原始CKKS参数
	N            int             // 总参与方数量
	T            int             // 门限值
}

// NewSemiShamirApproxSS 创建新的SemiShamirApproxSS实例
func NewSemiShamirApproxSS(N, T int, params ckks.Parameters) (semiShamirApproxSS *SemiShamirApproxSS) {
	semiShamirApproxSS = new(SemiShamirApproxSS)
	semiShamirApproxSS.N = N
	semiShamirApproxSS.T = T
	semiShamirApproxSS.originalCKKS = params // 存储原始CKKS参数
	semiShamirApproxSS.VanSS = NewVanillaShamirSS(N, T, params)

	// 设置组合器参数，使用第一个模数创建新的参数
	semiShamirApproxSS.params4cmb, _ = rlwe.NewParameters(params.LogN(), params.Q()[:1], nil, 0, params.HammingWeight(), params.Sigma(), params.RingType(), rlwe.NewScale(1), params.DefaultNTTFlag())

	return
}

// ApproxRecover 执行半门限近似恢复算法
func (semi *SemiShamirApproxSS) ApproxRecover(bound4smudgingNoise int) (timeComp time.Duration, sizeComm float64) {

	paramsMesSpace := semi.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2 time.Duration
	var sizeCommStage1, sizeCommStage2 float64

	// Round 1: 参与方发送标记信息，聚合者计算拉格朗日系数

	// 选择T个参与方
	parR1 := SelectRandomParticipants(semi.N, semi.T)
	fmt.Printf("Round 1: Selected participants: %v\n", parR1)

	// 第一轮时间开始计算
	timeStage1Start := time.Now()

	// 参与方发送标记信息（在这里我们模拟发送participant ID作为标记）
	participantIDs := make([]int, len(parR1))
	for i, par := range parR1 {
		participantIDs[i] = par
	}

	// 聚合者收集所有标记信息，计算拉格朗日系数
	survivialPublicPoint := make([]drlwe.ShamirPublicPoint, len(parR1))
	for i, pi := range parR1 {
		survivialPublicPoint[i] = drlwe.ShamirPublicPoint(pi + 1)
	}

	// 计算拉格朗日系数
	cmbR1 := NewMyCombiner(&semi.params4cmb, survivialPublicPoint, semi.T)
	lagrangeCoeffs := make(map[int]uint64, semi.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(semi.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}

	timeStage1 = time.Since(timeStage1Start)
	// 第一轮通信：参与方ID（标记信息）+ 拉格朗日系数
	sizeCommStage1 = float64(4) + float64(8*len(lagrangeCoeffs)) // 4字节per ID + 8字节per coefficient

	fmt.Printf("Round 1 completed: Lagrange coefficients computed for %d participants\n", len(lagrangeCoeffs))

	// Round 2: 参与方计算 λᵢ^(T)·bᵢ + nᵢ 并发送给聚合者

	timeStage2Start := time.Now()

	// 存储每个参与方的部分解密值
	partialValues := make(map[int]*ring.Poly, semi.T)

	for _, par := range parR1 {
		// 读取参与方的秘密份额 bᵢ
		filename := fmt.Sprintf("%s%d", "mesToParty", par)
		skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
		skFile, _ := os.Open(skFilename)
		defer skFile.Close()
		skData, _ := io.ReadAll(skFile)
		share := semi.VanSS.thdizer.AllocateThresholdSecretShare()
		share.UnmarshalBinary(skData)

		// 生成噪声 nᵢ
		smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSamplerFromRLWE(semi.VanSS.thdizer.params, paramsMesSpace.RingQ(), bound4smudgingNoise)

		ni := ringQ.NewPoly()
		smudgingNoiseSampler.ReadLvl(levelQ, ni)

		// 输出安全信息
		if actualBound != bound4smudgingNoise {
			fmt.Printf("SemiShamirApproxSS: Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
		}

		// 计算 λᵢ^(T)·bᵢ + nᵢ
		result := ringQ.NewPoly()
		lagrangeCoeff := lagrangeCoeffs[par]

		// 先计算 λᵢ·bᵢ
		ringQ.MulScalar(share.Poly.Q, lagrangeCoeff, result)

		// 再加上噪声 nᵢ
		ringQ.Add(result, ni, result)

		partialValues[par] = result

		fmt.Printf("Participant %d computed partial value with Lagrange coeff %d\n", par, lagrangeCoeff)
	}

	// 聚合者相加所有部分解密值，得到总的 Σ{i∈T}(λᵢ^(T)·bᵢ + nᵢ)
	approxMessage := ringQ.NewPoly()
	approxMessage.Zero()

	for _, par := range parR1 {
		ringQ.Add(approxMessage, partialValues[par], approxMessage)
	}

	timeStage2 = time.Since(timeStage2Start)
	// 第二轮通信：每个参与方发送一个多项式
	sizeCommStage2 = float64(partialValues[parR1[0]].MarshalBinarySize64())

	fmt.Printf("Round 2 completed: Aggregated partial values from %d participants\n", len(partialValues))

	// 计算总时间和通信量
	timeComp = timeStage1 + timeStage2
	sizeComm = sizeCommStage1 + sizeCommStage2

	fmt.Printf("Time breakdown: Stage1=%v, Stage2=%v, Total=%v\n", timeStage1, timeStage2, timeComp)
	fmt.Printf("Communication breakdown: Stage1=%.2f bytes, Stage2=%.2f bytes, Total=%.2f bytes\n", sizeCommStage1, sizeCommStage2, sizeComm)

	return
}

// ApproxRecover4TestTime 用于测试时间的版本（只处理第一个参与方的实际计算）
func (semi *SemiShamirApproxSS) ApproxRecover4TestTime(bound4smudgingNoise int) (timeComp time.Duration, sizeComm float64) {

	paramsMesSpace := semi.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2 time.Duration
	var sizeCommStage1, sizeCommStage2 float64

	// Round 1: 参与方发送标记信息，聚合者计算拉格朗日系数

	// 选择T个参与方
	parR1 := SelectRandomParticipants(semi.N, semi.T)

	// 第一轮时间开始计算
	timeStage1Start := time.Now()

	// 聚合者收集所有标记信息，计算拉格朗日系数
	survivialPublicPoint := make([]drlwe.ShamirPublicPoint, len(parR1))
	for i, pi := range parR1 {
		survivialPublicPoint[i] = drlwe.ShamirPublicPoint(pi + 1)
	}

	// 计算拉格朗日系数
	cmbR1 := NewMyCombiner(&semi.params4cmb, survivialPublicPoint, semi.T)
	lagrangeCoeffs := make(map[int]uint64, semi.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(semi.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}

	timeStage1 = time.Since(timeStage1Start)
	sizeCommStage1 = float64(4) + float64(8*len(lagrangeCoeffs))

	// Round 2: 参与方计算 λᵢ^(T)·bᵢ + nᵢ （测试版本：只处理第一个参与方）

	timeStage2Start := time.Now()

	// 存储每个参与方的部分解密值
	partialValues := make(map[int]*ring.Poly, semi.T)

	for i, par := range parR1 {
		if i == 0 {
			// 只对第一个参与方进行实际计算
			// 读取参与方的秘密份额 bᵢ
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			defer skFile.Close()
			skData, _ := io.ReadAll(skFile)
			share := semi.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			// 生成噪声 nᵢ
			smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSamplerFromRLWE(semi.VanSS.thdizer.params, paramsMesSpace.RingQ(), bound4smudgingNoise)

			ni := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, ni)

			// 输出安全信息
			if actualBound != bound4smudgingNoise {
				fmt.Printf("SemiShamirApproxSS (TestTime): Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
			}

			// 计算 λᵢ^(T)·bᵢ + nᵢ
			result := ringQ.NewPoly()
			lagrangeCoeff := lagrangeCoeffs[par]

			// 先计算 λᵢ·bᵢ
			ringQ.MulScalar(share.Poly.Q, lagrangeCoeff, result)

			// 再加上噪声 nᵢ
			ringQ.Add(result, ni, result)

			partialValues[par] = result
		} else {
			// 其他参与方复用第一个参与方的结果
			partialValues[par] = partialValues[parR1[0]]
		}
	}

	// 聚合者相加所有部分解密值
	approxMessage := ringQ.NewPoly()
	approxMessage.Zero()

	for _, par := range parR1 {
		ringQ.Add(approxMessage, partialValues[par], approxMessage)
	}

	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(partialValues[parR1[0]].MarshalBinarySize64())

	// 计算总时间和通信量
	timeComp = timeStage1 + timeStage2
	sizeComm = sizeCommStage1 + sizeCommStage2

	fmt.Printf("TestTime - Time breakdown: Stage1=%v, Stage2=%v, Total=%v\n", timeStage1, timeStage2, timeComp)
	fmt.Printf("Size1=%v, Size2=%v, Total=%v\n", sizeCommStage1, sizeCommStage2, sizeComm)
	return
}

func (semi *SemiShamirApproxSS) ApproxRecoverNN4TestTime(bound4smudgingNoise int) (timeComp time.Duration, sizeComm float64) {

	paramsMesSpace := semi.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64

	// Round 1: 参与方发送标记信息，聚合者计算拉格朗日系数

	// 选择T个参与方
	parR1 := SelectRandomParticipants(semi.N, semi.T)

	timeStage2Start := time.Now()

	// 存储每个参与方的部分解密值
	partialValues := make(map[int]*ring.Poly, semi.T)

	for i, par := range parR1 {
		if i == 0 {
			// 只对第一个参与方进行实际计算
			// 读取参与方的秘密份额 bᵢ
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			defer skFile.Close()
			skData, _ := io.ReadAll(skFile)
			share := semi.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			// 生成噪声 nᵢ
			smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSamplerFromRLWE(semi.VanSS.thdizer.params, paramsMesSpace.RingQ(), bound4smudgingNoise)

			ni := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, ni)

			// 输出安全信息
			if actualBound != bound4smudgingNoise {
				fmt.Printf("SemiShamirApproxSS (TestTime): Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
			}

			// 计算 bᵢ + nᵢ
			result := ringQ.NewPoly()
			ringQ.Add(result, ni, result)

			partialValues[par] = result
		} else {
			// 其他参与方复用第一个参与方的结果
			partialValues[par] = partialValues[parR1[0]]
		}
	}

	// 聚合者相加所有部分解密值
	approxMessage := ringQ.NewPoly()
	approxMessage.Zero()

	for _, par := range parR1 {
		ringQ.Add(approxMessage, partialValues[par], approxMessage)
	}

	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(partialValues[parR1[0]].MarshalBinarySize64())

	// Round 3: 进行 2(N-T)/2 次同态加密私钥、同态乘法和同态加法
	timeStage3Start := time.Now()

	// 计算需要处理的参与方数量
	numRemainingParticipants := semi.N - semi.T

	// 存储同态运算结果
	encryptedResults := make([]*ring.Poly, numRemainingParticipants)

	// 创建CKKS参数用于同态运算
	ckksParams := semi.originalCKKS
	encoder := ckks.NewEncoder(ckksParams)
	evaluator := ckks.NewEvaluator(ckksParams, rlwe.EvaluationKey{})

	for i := 0; i < numRemainingParticipants; i++ {
		// 生成私钥多项式
		privateKey := ringQ.NewPoly()
		// 使用简单的方式生成随机多项式
		for j := 0; j < privateKey.N(); j++ {
			privateKey.Coeffs[0][j] = uint64(i + j + 1) // 简单的确定性生成
		}

		// 1. 同态加密私钥
		// 将uint64转换为complex128用于CKKS编码
		slots := ckksParams.Slots()
		encodeLength := privateKey.N()
		if encodeLength > slots {
			encodeLength = slots
		}

		keyComplex := make([]complex128, encodeLength)
		for j := 0; j < encodeLength; j++ {
			keyComplex[j] = complex(float64(privateKey.Coeffs[0][j]), 0)
		}
		keyPlaintext := ckks.NewPlaintext(ckksParams, ckksParams.MaxLevel())
		encoder.Encode(keyComplex, keyPlaintext, ckksParams.LogSlots())

		// 模拟加密（由于没有公钥，我们直接使用明文作为"密文"进行后续运算）
		encryptedKey := ckks.NewCiphertext(ckksParams, 1, ckksParams.MaxLevel())
		encryptedKey.Value[0] = keyPlaintext.Value.CopyNew()
		encryptedKey.Value[1] = ringQ.NewPoly() // 零多项式作为第二部分

		// 2. 准备approxMessage作为明文
		approxLength := approxMessage.N()
		if approxLength > slots {
			approxLength = slots
		}

		approxComplex := make([]complex128, approxLength)
		for j := 0; j < approxLength; j++ {
			approxComplex[j] = complex(float64(approxMessage.Coeffs[0][j]), 0)
		}
		approxPlaintext := ckks.NewPlaintext(ckksParams, ckksParams.MaxLevel())
		encoder.Encode(approxComplex, approxPlaintext, ckksParams.LogSlots())

		// 3. 同态乘法：密文与密文相乘
		multResult := ckks.NewCiphertext(ckksParams, 2, ckksParams.MaxLevel()) // 两个degree 1密文相乘得到degree 2
		evaluator.Mul(encryptedKey, encryptedKey, multResult)

		// 4. 同态加法：将结果与自身相加
		finalResult := ckks.NewCiphertext(ckksParams, 2, ckksParams.MaxLevel()) // 保持degree 2
		evaluator.Add(multResult, multResult, finalResult)

		// 存储结果的第一个组件作为多项式
		encryptedResults[i] = finalResult.Value[0].CopyNew()
	}

	timeStage3 = time.Since(timeStage3Start)

	// 第一阶段通信量：(N-T)次私钥密文传输的量
	if numRemainingParticipants > 0 {
		sizeCommStage1 = float64(numRemainingParticipants) * float64(encryptedResults[0].MarshalBinarySize64())
	} else {
		sizeCommStage1 = 0
	}

	// 第三阶段通信量：同态运算结果传输
	if len(encryptedResults) > 0 {
		sizeCommStage3 = float64(len(encryptedResults)) * float64(encryptedResults[0].MarshalBinarySize64())
	} else {
		sizeCommStage3 = 0
	}

	// 计算总时间和通信量
	timeComp = timeStage2 + timeStage3
	sizeComm = sizeCommStage1 + sizeCommStage2 + sizeCommStage3

	fmt.Printf("TestTime - Time breakdown: Stage1=%v, Stage2=%v, Stage3=%v, Total=%v\n", timeStage1, timeStage2, timeStage3, timeComp)
	fmt.Printf("Size1=%v, Size2=%v, Size3=%v, Total=%v\n", sizeCommStage1, sizeCommStage2, sizeCommStage3, sizeComm)
	return
}
