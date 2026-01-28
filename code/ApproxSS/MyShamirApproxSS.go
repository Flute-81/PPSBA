package ApproxSS

import (
	"fmt"
	"io"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type MyShamirApproxSS struct {
	VanSS *VanillaShamirSS

	params4DoubleEncryption ckks.Parameters
	params4cmb              rlwe.Parameters
	f                       *skEncryptionCKKS
	originalParams          ckks.Parameters // 存储原始CKKS参数
	ek1All                  []*rlwe.SecretKey
	ek2All                  []*rlwe.SecretKey
	N                       int
	T                       int
}

func NewMyShamirApproxSS(N, T int, params ckks.Parameters) (myShamirApproxSS *MyShamirApproxSS) {
	myShamirApproxSS = new(MyShamirApproxSS)
	myShamirApproxSS.N = N
	myShamirApproxSS.T = T
	myShamirApproxSS.originalParams = params // 保存原始CKKS参数
	myShamirApproxSS.VanSS = NewVanillaShamirSS(N, T, params)
	// 对于当前实现，我们直接复用原始 CKKS 参数进行“双重加密”部分，
	// 避免额外的参数放大和素数搜索带来的时间开销。

	myShamirApproxSS.params4DoubleEncryption = findLargerParameters(params, N)

	myShamirApproxSS.f = NewSKencryption(myShamirApproxSS.params4DoubleEncryption)
	myShamirApproxSS.params4cmb, _ = rlwe.NewParameters(params.LogN(), params.Q()[:1], nil, 0, params.HammingWeight(), params.Sigma(), params.RingType(), rlwe.NewScale(1), params.DefaultNTTFlag())

	myShamirApproxSS.ek1All = make([]*rlwe.SecretKey, N)
	myShamirApproxSS.ek2All = make([]*rlwe.SecretKey, N)

	//Still need to generate and share the encryption key
	return
}

func findLargerParameters(params ckks.Parameters, N int) (params4DoubleEncryption ckks.Parameters) {

	return params
}

// func (myShamirApproxSS *MyShamirApproxSS) Preprocessing(bound4smudgingNoise int) (timeComp time.Duration, sizeComm float64) {

// 	paramsMesSpace := myShamirApproxSS.VanSS.thdizer.params
// 	ringQ := paramsMesSpace.RingQ()
// 	levelQ := 0

// 	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
// 	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64

// 	//Round 1 Starts

// 	//Samples the polynomial a in our paper
// 	CRS, _ := utils.NewPRNG()
// 	a := TagGen(myShamirApproxSS.params4DoubleEncryption, CRS)

// 	//Parties run PartyR1
// 	for i := 0; i < myShamirApproxSS.VanSS.N; i++ {
// 		//Read the share from file
// 		filename := fmt.Sprintf("%s%d", "mesToParty", i)
// 		skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
// 		skFile, _ := os.Open(skFilename)
// 		defer skFile.Close()
// 		skData, _ := io.ReadAll(skFile)
// 		share := myShamirApproxSS.VanSS.thdizer.AllocateThresholdSecretShare()
// 		share.UnmarshalBinary(skData)

// 		timeStage1Start := time.Now()
// 		prng_i, _ := utils.NewPRNG()
// 		smudgingNoiseSampler := ring.NewGaussianSampler(prng_i, paramsMesSpace.RingQ(), float64(bound4smudgingNoise)/6, bound4smudgingNoise)

// 		ni := ringQ.NewPoly()
// 		smudgingNoiseSampler.ReadLvl(levelQ, ni)
// 		CTni_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)

// 		CTni_all[i] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek2All[i], ni)
// 		timeStage1 = time.Since(timeStage1Start)
// 		sizeCommStage1 = float64(CTni_all[i].Value[0].MarshalBinarySize64())
// 	}

// }

func (myShamirApproxSS *MyShamirApproxSS) ApproxRecover(bound4smudgingNoise int) (isSuc bool, timeCompOnce time.Duration, timeCompNotOnce time.Duration, sizeCommOnce float64, sizeCommNotOnce float64) {

	paramsMesSpace := myShamirApproxSS.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
	var timeNotOnce1, timeNotOnce2 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64
	var sizeCommNotOnce1 float64

	//Round 1 Starts

	//Samples the polynomial a in our paper
	CRS, _ := utils.NewPRNG()
	a := TagGen(myShamirApproxSS.params4DoubleEncryption, CRS)
	//Select the participants in round 1
	parR1 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	//The output of parties and aggregator in Round 1
	CTni_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	CTsi_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	lagrangeCoeffs := make(map[int]uint64, myShamirApproxSS.VanSS.T)

	//Parties run PartyR1
	for _, par := range parR1 {
		//Read the share to be encrypted from file
		filename := fmt.Sprintf("%s%d", "mesToParty", par)
		skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
		skFile, _ := os.Open(skFilename)
		defer skFile.Close()
		skData, _ := io.ReadAll(skFile)
		share := myShamirApproxSS.VanSS.thdizer.AllocateThresholdSecretShare()
		share.UnmarshalBinary(skData)

		myShamirApproxSS.ek1All[par], timeStage1 = myShamirApproxSS.f.GenerateEKThenShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "s")
		myShamirApproxSS.ek2All[par], _ = myShamirApproxSS.f.GenerateEKThenShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "n")
		timeStage1 = timeStage1 * 2
		sizeCommStage1 = float64(myShamirApproxSS.ek1All[par].Value.MarshalBinarySize64()*myShamirApproxSS.N) * 2

		timeNotOnce1Start := time.Now()
		// 使用安全的smudging noise采样器
		smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSampler(myShamirApproxSS.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)

		ni := ringQ.NewPoly()
		smudgingNoiseSampler.ReadLvl(levelQ, ni)

		// 输出安全信息
		if actualBound != bound4smudgingNoise {
			fmt.Printf("MyShamirApproxSS: Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
		}

		CTsi_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek1All[par], share.Poly.Q)
		CTni_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek2All[par], ni)
		timeNotOnce1 = time.Since(timeNotOnce1Start)
		sizeCommNotOnce1 = float64(CTsi_all[par].Value[0].MarshalBinarySize64()) * 2
	}

	//AggregatorR1
	timeStage2Start := time.Now()
	survivialPublicPoint := []drlwe.ShamirPublicPoint{}
	//Collect survivialClient
	for _, pi := range parR1 {
		survivialPublicPoint = append(survivialPublicPoint, drlwe.ShamirPublicPoint(pi+1))
	}
	cmbR1 := NewMyCombiner(&myShamirApproxSS.params4cmb, survivialPublicPoint, myShamirApproxSS.VanSS.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(myShamirApproxSS.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}
	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(8 * len(lagrangeCoeffs))
	//Round 1 Ends!

	//Round 2 Starts!

	//Select the participants in round 1
	parR2 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	//The output of parties in Round 1
	dkShare_All := make(map[int]drlwe.ShamirSecretShare, myShamirApproxSS.VanSS.T)

	//Parties run PartyR2
	for _, par := range parR2 {

		dkShare_All[par], timeStage3 = myShamirApproxSS.f.GenerateDKShareFile(parR1, lagrangeCoeffs, par)

		sizeCommStage3 = float64(dkShare_All[par].Poly.Q.MarshalBinarySize64())
	}

	//AggregatorR2
	timeStage4Start := time.Now()
	DK := myShamirApproxSS.f.GenerateDK(parR2, dkShare_All)
	timeStage4 = time.Since(timeStage4Start)

	timeNotOnce2Start := time.Now()
	approxMessage := ringQ.NewPoly()
	myShamirApproxSS.f.FEDecFinal(DK, CTni_all, CTsi_all, lagrangeCoeffs, parR1, approxMessage)
	timeNotOnce2 = time.Since(timeNotOnce2Start)

	//Round 2 Ends!

	timeCompOnce = timeStage1 + timeStage2 + timeStage3 + timeStage4
	fmt.Println(timeStage1, timeCompOnce)
	timeCompNotOnce = timeNotOnce1 + timeNotOnce2
	sizeCommOnce = sizeCommStage1 + sizeCommStage2 + sizeCommStage3
	sizeCommNotOnce = sizeCommNotOnce1
	isSuc = TestResult(approxMessage, myShamirApproxSS.VanSS.secret.Value.Q, uint64(myShamirApproxSS.VanSS.T*bound4smudgingNoise), myShamirApproxSS.VanSS.thdizer.params.Q()[0:1])
	return

}

func (myShamirApproxSS *MyShamirApproxSS) ApproxRecover4TestTime(bound4smudgingNoise int) (timeCompOnce time.Duration, timeCompNotOnce time.Duration, sizeCommOnce float64, sizeCommNotOnce float64) {

	paramsMesSpace := myShamirApproxSS.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
	var timeNotOnce1, timeNotOnce2 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64
	var sizeCommNotOnce1 float64

	//Round 1 Starts

	//Samples the polynomial a in our paper
	CRS, _ := utils.NewPRNG()
	a := TagGen(myShamirApproxSS.params4DoubleEncryption, CRS)
	//Select the participants in round 1
	parR1 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	//The output of parties and aggregator in Round 1
	CTni_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	CTsi_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	lagrangeCoeffs := make(map[int]uint64, myShamirApproxSS.VanSS.T)

	//Parties run PartyR1
	for i, par := range parR1 {
		if i == 0 {
			//Read the share to be encrypted from file
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			defer skFile.Close()
			skData, _ := io.ReadAll(skFile)
			share := myShamirApproxSS.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			myShamirApproxSS.ek1All[par], timeStage1 = myShamirApproxSS.f.GenerateEKThenShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "s")
			myShamirApproxSS.ek2All[par], _ = myShamirApproxSS.f.GenerateEKThenShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "n")
			timeStage1 = timeStage1 * 2
			sizeCommStage1 = float64(myShamirApproxSS.ek1All[par].Value.MarshalBinarySize64()*myShamirApproxSS.N) * 2

			timeNotOnce1Start := time.Now()
			// 使用安全的smudging noise采样器
			smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSampler(myShamirApproxSS.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)

			ni := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, ni)

			// 输出安全信息
			if actualBound != bound4smudgingNoise {
				fmt.Printf("MyShamirApproxSS (TestTime): Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
			}

			CTsi_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek1All[par], share.Poly.Q)
			CTni_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek2All[par], ni)
			timeNotOnce1 = time.Since(timeNotOnce1Start)
			sizeCommNotOnce1 = float64(CTsi_all[par].Value[0].MarshalBinarySize64()) * 2
		} else {
			myShamirApproxSS.ek1All[par] = myShamirApproxSS.ek1All[parR1[0]]
			myShamirApproxSS.ek2All[par] = myShamirApproxSS.ek2All[parR1[0]]
			CTsi_all[par] = CTsi_all[parR1[0]]
			CTni_all[par] = CTni_all[parR1[0]]
		}
	}

	//AggregatorR1
	timeStage2Start := time.Now()
	survivialPublicPoint := []drlwe.ShamirPublicPoint{}
	//Collect survivialClient
	for _, pi := range parR1 {
		survivialPublicPoint = append(survivialPublicPoint, drlwe.ShamirPublicPoint(pi+1))
	}
	cmbR1 := NewMyCombiner(&myShamirApproxSS.params4cmb, survivialPublicPoint, myShamirApproxSS.VanSS.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(myShamirApproxSS.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}
	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(8 * len(lagrangeCoeffs))
	//Round 1 Ends!

	//Round 2 Starts!

	//Select the participants in round 2
	parR2 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	//The output of parties in Round 2
	dkShare_All := make(map[int]drlwe.ShamirSecretShare, myShamirApproxSS.VanSS.T)

	//Parties run PartyR2
	for i, par := range parR2 {
		if i == 0 {
			dkShare_All[par], timeStage3 = myShamirApproxSS.f.GenerateDKShareFile4TestTime(parR1, lagrangeCoeffs, par)
			sizeCommStage3 = float64(dkShare_All[par].Poly.Q.MarshalBinarySize64())
		} else {
			dkShare_All[par] = dkShare_All[parR2[0]]
		}

	}

	//AggregatorR2
	timeStage4Start := time.Now()
	DK := myShamirApproxSS.f.GenerateDK(parR2, dkShare_All)
	timeStage4 = time.Since(timeStage4Start)

	timeNotOnce2Start := time.Now()
	approxMessage := ringQ.NewPoly()
	myShamirApproxSS.f.FEDecFinal(DK, CTni_all, CTsi_all, lagrangeCoeffs, parR1, approxMessage)
	timeNotOnce2 = time.Since(timeNotOnce2Start)

	//Round 2 Ends!

	timeCompOnce = timeStage1 + timeStage2 + timeStage3 + timeStage4
	fmt.Println(timeStage1, timeCompOnce)
	timeCompNotOnce = timeNotOnce1 + timeNotOnce2
	sizeCommOnce = sizeCommStage1 + sizeCommStage2 + sizeCommStage3
	sizeCommNotOnce = sizeCommNotOnce1

	fmt.Println("timeStageMyShamir：", timeStage1, timeStage2, timeStage3, timeStage4)
	fmt.Println("timeNotOnceMyShamir： ", timeNotOnce1, timeNotOnce2)
	fmt.Println("sizeCommOnceMyShamir： ", sizeCommStage1, sizeCommStage2, sizeCommStage3)
	fmt.Println("sizeNotOnceMyShamir： ", sizeCommNotOnce1)
	return

}

func (myShamirApproxSS *MyShamirApproxSS) ApproxSecRecover4TestTime(bound4smudgingNoise int) (timeCompOnce time.Duration, timeCompNotOnce time.Duration, sizeCommOnce float64, sizeCommNotOnce float64) {

	paramsMesSpace := myShamirApproxSS.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
	var timeNotOnce1, timeNotOnce2 time.Duration
	var sizeCommStage1, sizeCommStage2, sizeCommStage3 float64
	var sizeCommNotOnce1 float64

	//Round 1 Starts

	//Samples the polynomial a in our paper
	CRS, _ := utils.NewPRNG()
	a := TagGen(myShamirApproxSS.params4DoubleEncryption, CRS)
	//Select the participants in round 1
	parR1 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	//The output of parties and aggregator in Round 1
	CTni_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	CTsi_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	lagrangeCoeffs := make(map[int]uint64, myShamirApproxSS.VanSS.T)

	//Parties run PartyR1
	for i, par := range parR1 {
		if i == 0 {
			//Read the share to be encrypted from file
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			defer skFile.Close()
			skData, _ := io.ReadAll(skFile)
			share := myShamirApproxSS.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			// 使用新的函数：生成私钥份额，转换为a·ek+e'形式，然后分发
			myShamirApproxSS.ek1All[par], timeStage1 = myShamirApproxSS.f.GenerateEKThenConvertAndShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "s", a, bound4smudgingNoise)
			myShamirApproxSS.ek2All[par], _ = myShamirApproxSS.f.GenerateEKThenConvertAndShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "n", a, bound4smudgingNoise)
			timeStage1 = timeStage1 * 2

			// 计算转换后份额的大小：a·ek1+e1' 和 a·ek2+e2'
			// 生成噪声来计算转换后份额的大小
			smudgingNoiseSampler, _ := CreateSecureSmudgingNoiseSampler(myShamirApproxSS.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)
			e1Prime := ringQ.NewPoly()
			e2Prime := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, e1Prime)
			smudgingNoiseSampler.ReadLvl(levelQ, e2Prime)

			// 计算转换后的份额
			convertedShare1 := ringQ.NewPoly()
			convertedShare2 := ringQ.NewPoly()
			ringQ.MulCoeffsMontgomeryLvl(levelQ, a, myShamirApproxSS.ek1All[par].Value.Q, convertedShare1)
			ringQ.AddLvl(levelQ, convertedShare1, e1Prime, convertedShare1)
			ringQ.MulCoeffsMontgomeryLvl(levelQ, a, myShamirApproxSS.ek2All[par].Value.Q, convertedShare2)
			ringQ.AddLvl(levelQ, convertedShare2, e2Prime, convertedShare2)

			// 记录转换后份额的大小
			sizeCommStage1 = float64(convertedShare1.MarshalBinarySize64()+convertedShare2.MarshalBinarySize64()) * float64(myShamirApproxSS.N)

			timeNotOnce1Start := time.Now()
			// 使用安全的smudging noise采样器
			smudgingNoiseSampler, actualBound := CreateSecureSmudgingNoiseSampler(myShamirApproxSS.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)

			ni := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, ni)

			// 输出安全信息
			if actualBound != bound4smudgingNoise {
				fmt.Printf("MyShamirApproxSS (TestTime): Using secure bound %d (original: %d)\n", actualBound, bound4smudgingNoise)
			}

			CTsi_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek1All[par], share.Poly.Q)
			CTni_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek2All[par], ni)
			timeNotOnce1 = time.Since(timeNotOnce1Start)
			sizeCommNotOnce1 = float64(CTsi_all[par].Value[0].MarshalBinarySize64()) * 2
		} else {
			myShamirApproxSS.ek1All[par] = myShamirApproxSS.ek1All[parR1[0]]
			myShamirApproxSS.ek2All[par] = myShamirApproxSS.ek2All[parR1[0]]
			CTsi_all[par] = CTsi_all[parR1[0]]
			CTni_all[par] = CTni_all[parR1[0]]
		}
	}

	//AggregatorR1
	timeStage2Start := time.Now()
	survivialPublicPoint := []drlwe.ShamirPublicPoint{}
	//Collect survivialClient
	for _, pi := range parR1 {
		survivialPublicPoint = append(survivialPublicPoint, drlwe.ShamirPublicPoint(pi+1))
	}
	cmbR1 := NewMyCombiner(&myShamirApproxSS.params4cmb, survivialPublicPoint, myShamirApproxSS.VanSS.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(myShamirApproxSS.params4cmb.RingQP(), cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}
	timeStage2 = time.Since(timeStage2Start)
	sizeCommStage2 = float64(8 * len(lagrangeCoeffs))
	//Round 1 Ends!

	//Round 2 Starts!

	//Select the participants in round 2
	parR2 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	//The output of parties in Round 2
	dkShare_All := make(map[int]drlwe.ShamirSecretShare, myShamirApproxSS.VanSS.T)

	//Parties run PartyR2
	for i, par := range parR2 {
		if i == 0 {
			dkShare_All[par], timeStage3 = myShamirApproxSS.f.GenerateDKShareFile4TestTime(parR1, lagrangeCoeffs, par)
			// 修改：计算新的重构份额 λ_i·(a·ek_{i,1}+e_{i,1}') + (a·ek_{i,2}+e_{i,2}')
			// 生成新的噪声e1'和e2'（与阶段1中相同的采样器）
			smudgingNoiseSampler, _ := CreateSecureSmudgingNoiseSampler(myShamirApproxSS.originalParams, paramsMesSpace.RingQ(), bound4smudgingNoise)
			e1Prime := ringQ.NewPoly()
			e2Prime := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, e1Prime)
			smudgingNoiseSampler.ReadLvl(levelQ, e2Prime)

			// 计算转换后的份额：a·ek1 + e1' 和 a·ek2 + e2'
			convertedShare1 := ringQ.NewPoly()
			convertedShare2 := ringQ.NewPoly()

			// a·ek1 + e1'
			ringQ.MulCoeffsMontgomeryLvl(levelQ, a, myShamirApproxSS.ek1All[parR1[0]].Value.Q, convertedShare1)
			ringQ.AddLvl(levelQ, convertedShare1, e1Prime, convertedShare1)

			// a·ek2 + e2'
			ringQ.MulCoeffsMontgomeryLvl(levelQ, a, myShamirApproxSS.ek2All[parR1[0]].Value.Q, convertedShare2)
			ringQ.AddLvl(levelQ, convertedShare2, e2Prime, convertedShare2)

			// 计算新的重构份额：λ_i·(a·ek_{i,1}+e_{i,1}') + (a·ek_{i,2}+e_{i,2}')
			newReconShare := ringQ.NewPoly()
			lambdaScalar := ringQ.NewRNSScalarFromUInt64(lagrangeCoeffs[parR1[0]])

			// λ_i·(a·ek_{i,1}+e_{i,1}')
			ringQ.MulRNSScalarMontgomery(convertedShare1, lambdaScalar, newReconShare)

			// 加上 (a·ek_{i,2}+e_{i,2}')
			ringQ.AddLvl(levelQ, newReconShare, convertedShare2, newReconShare)

			// 修改：记录新重构份额的大小
			sizeCommStage3 = float64(newReconShare.MarshalBinarySize64())
		} else {
			dkShare_All[par] = dkShare_All[parR2[0]]
		}

	}

	//AggregatorR2
	timeStage4Start := time.Now()
	DK := myShamirApproxSS.f.GenerateDK(parR2, dkShare_All)
	timeStage4 = time.Since(timeStage4Start)

	timeNotOnce2Start := time.Now()
	approxMessage := ringQ.NewPoly()
	myShamirApproxSS.f.FEDecFinal(DK, CTni_all, CTsi_all, lagrangeCoeffs, parR1, approxMessage)
	timeNotOnce2 = time.Since(timeNotOnce2Start)

	//Round 2 Ends!

	timeCompOnce = timeStage1 + timeStage2 + timeStage3 + timeStage4
	fmt.Println(timeStage1, timeCompOnce)
	timeCompNotOnce = timeNotOnce1 + timeNotOnce2
	sizeCommOnce = sizeCommStage1 + sizeCommStage2 + sizeCommStage3
	sizeCommNotOnce = sizeCommNotOnce1

	fmt.Println("timeStageMyShamirSec：", timeStage1, timeStage2, timeStage3, timeStage4)
	fmt.Println("timeNotOnceMyShamirSec： ", timeNotOnce1, timeNotOnce2)
	return

}
