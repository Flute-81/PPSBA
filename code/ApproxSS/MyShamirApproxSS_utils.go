package ApproxSS

import (
	"fmt"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func TagGen(params ckks.Parameters, prng *utils.KeyedPRNG) *ring.Poly {
	tag := params.RingQP().NewPoly()
	uniformSampler := ringqp.NewUniformSampler(prng, *params.RingQP())
	uniformSampler.ReadLvl(params.MaxLevel(), -1, tag)
	return tag.Q
}

type skEncryptionCKKS struct {
	ckksParams          ckks.Parameters
	gaussianSamplerQ4FE *ring.GaussianSampler
	uniformSamplerQ4FE  ringqp.UniformSampler
	thr4FE              *MyThresholdizer
	encoder             ckks.Encoder
}

// func NewFE(params rlwe.Parameters, T uint64) *skEncryptionCKKS {
// 	f := new(FE)
// 	f.params = params
// 	f.T = T
// 	prng, _ := utils.NewPRNG()
// 	f.gaussianSamplerQ4FE = ring.NewGaussianSampler(prng, f.params.RingQ(), f.params.Sigma(), int(6*f.params.Sigma()))
// 	f.uniformSamplerQ4FE = ring.NewUniformSampler(prng, f.params.RingQ())
// 	f.thr4FE = drlwe.NewThresholdizer(f.params)
// 	deltaInt := f.params.Q()[0] / f.T
// 	*f.delta = f.params.RingQ().NewRNSScalarFromUInt64(deltaInt)
// 	f.scaler = bfv.NewRNSScaler(f.params.RingQ(), f.T)

// 	return f
// }

func NewSKencryption(params ckks.Parameters) (f *skEncryptionCKKS) {

	f = new(skEncryptionCKKS)
	f.ckksParams = params

	prng, _ := utils.NewPRNG()
	f.gaussianSamplerQ4FE = ring.NewGaussianSampler(prng, f.ckksParams.RingQ(), f.ckksParams.Sigma(), int(6*f.ckksParams.Sigma()))
	f.uniformSamplerQ4FE = ringqp.NewUniformSampler(prng, *f.ckksParams.RingQP())
	f.thr4FE = NewMyThresholdizer(f.ckksParams.Parameters)
	f.encoder = ckks.NewEncoder(f.ckksParams)
	return
}

func (f *skEncryptionCKKS) EncPolyCoeff(tag *ring.Poly, sk *rlwe.SecretKey, mes *ring.Poly) *rlwe.Ciphertext {

	mesVector := make([]uint64, mes.N())
	for i := 0; i < mes.N(); i++ {
		mesVector[i] = mes.Buff[i]
	}

	ct := f.FEEnc(tag, sk, mesVector)
	return ct
}

func (f *skEncryptionCKKS) FEEnc(tag *ring.Poly, sk *rlwe.SecretKey, mesVector []uint64) *rlwe.Ciphertext {

	pt := ckks.NewPlaintext(f.ckksParams, f.ckksParams.MaxLevel())
	// 将uint64转换为complex128用于CKKS编码
	// CKKS只能编码slots个值，不能超过
	slots := f.ckksParams.Slots()
	encodeLength := len(mesVector)
	if encodeLength > slots {
		encodeLength = slots
	}

	mesComplex := make([]complex128, encodeLength)
	for i := 0; i < encodeLength; i++ {
		mesComplex[i] = complex(float64(mesVector[i]), 0)
	}
	f.encoder.Encode(mesComplex, pt, f.ckksParams.LogSlots())

	ct := ckks.NewCiphertext(f.ckksParams, 1, pt.Level())
	ct.Value[1] = tag.CopyNew()
	c1 := ct.Value[1]

	ringQ := f.ckksParams.RingQ()
	levelQ := ct.Level()

	c0 := ct.Value[0]

	ringQ.MulCoeffsMontgomeryLvl(levelQ, c1, sk.Value.Q, c0) // c0 = NTT(sc1)
	ringQ.NegLvl(levelQ, c0, c0)                             // c0 = NTT(-sc1)

	buff := f.ckksParams.RingQ().NewPoly()
	if ct.IsNTT {
		f.gaussianSamplerQ4FE.ReadLvl(levelQ, buff) // e

		ringQ.NTTLvl(levelQ, buff, buff)   // NTT(e)
		ringQ.AddLvl(levelQ, c0, buff, c0) // c0 = NTT(-sc1 + e)
	} else {
		ringQ.InvNTTLvl(levelQ, c0, c0) // c0 = -sc1
		if ct.Degree() == 1 {
			ringQ.InvNTTLvl(levelQ, c1, c1) // c1 = c1
		}
		f.gaussianSamplerQ4FE.ReadAndAddLvl(levelQ, c0) // c0 = -sc1 + e
	}

	f.ckksParams.RingQ().AddLvl(ct.Level(), ct.Value[0], pt.Value, ct.Value[0])

	return ct

}

func (f *skEncryptionCKKS) GenerateEKfromShareFile(N, generator int, NorS string) (sk *rlwe.SecretKey) {

	sk = rlwe.NewSecretKey(f.ckksParams.Parameters)
	levelP := sk.LevelP()
	levelQ := sk.LevelQ()
	uShare := f.thr4FE.AllocateThresholdSecretShare()
	for j := 0; j < N; j++ {
		uShare.UnmarshalBinary(readSK4FEShare(generator, j, NorS))
		f.ckksParams.Parameters.RingQP().AddLvl(levelQ, levelP, sk.Value, uShare.Poly, sk.Value)
	}

	return sk
}

func (f *skEncryptionCKKS) GenerateEKThenShareToFile(N, T, generator int, NorS string) (*rlwe.SecretKey, time.Duration) {

	start := time.Now()
	keyGen4FE := ckks.NewKeyGenerator(f.ckksParams)
	sk := keyGen4FE.GenSecretKey()
	points := make([]drlwe.ShamirPublicPoint, N)
	for i := 0; i < N; i++ {
		points[i] = drlwe.ShamirPublicPoint(i + 1)
	}

	sk4FE1Share := f.thr4FE.GenShamirSecretShares(T, points, sk)

	timeTrue := time.Since(start)

	for i := 0; i < N; i++ {
		writeSK4FEShare(i, generator, NorS, *sk4FE1Share[i])
	}

	return sk, timeTrue
}

// GenerateEKThenConvertAndShareToFile 生成私钥份额，转换为a·ek+e'形式，然后分发
func (f *skEncryptionCKKS) GenerateEKThenConvertAndShareToFile(N, T, generator int, NorS string, a *ring.Poly, bound4smudgingNoise int) (*rlwe.SecretKey, time.Duration) {

	start := time.Now()
	keyGen4FE := ckks.NewKeyGenerator(f.ckksParams)
	sk := keyGen4FE.GenSecretKey()
	points := make([]drlwe.ShamirPublicPoint, N)
	for i := 0; i < N; i++ {
		points[i] = drlwe.ShamirPublicPoint(i + 1)
	}

	sk4FE1Share := f.thr4FE.GenShamirSecretShares(T, points, sk)

	// 生成新的噪声e'
	ringQ := f.ckksParams.RingQ()
	levelQ := 0
	smudgingNoiseSampler, _ := CreateSecureSmudgingNoiseSampler(f.ckksParams, ringQ, bound4smudgingNoise)

	ePrime := ringQ.NewPoly()
	smudgingNoiseSampler.ReadLvl(levelQ, ePrime)

	// 计算转换后的份额：a·ek + e'
	convertedShare := ringQ.NewPoly()
	ringQ.MulCoeffsMontgomeryLvl(levelQ, a, sk.Value.Q, convertedShare)
	ringQ.AddLvl(levelQ, convertedShare, ePrime, convertedShare)

	// 将转换后的份额写入文件
	for i := 0; i < N; i++ {
		// 为每个参与者计算转换后的份额
		participantConvertedShare := ringQ.NewPoly()
		ringQ.MulCoeffsMontgomeryLvl(levelQ, a, sk4FE1Share[i].Poly.Q, participantConvertedShare)
		ringQ.AddLvl(levelQ, participantConvertedShare, ePrime, participantConvertedShare)

		// 创建转换后的份额结构
		convertedShareStruct := f.thr4FE.AllocateThresholdSecretShare()
		convertedShareStruct.Poly.Q = participantConvertedShare
		convertedShareStruct.Poly.P = sk4FE1Share[i].Poly.P // 保持P部分不变

		writeSK4FEShare(i, generator, NorS, *convertedShareStruct)
	}

	timeTrue := time.Since(start)

	return sk, timeTrue
}

func (f *skEncryptionCKKS) GenerateDKShareNew(N int, eknShare []*drlwe.ShamirSecretShare, eks *drlwe.ShamirSecretShare) drlwe.ShamirSecretShare {

	dkShare := f.thr4FE.AllocateThresholdSecretShare()

	dkShare.Poly = f.ckksParams.RingQP().NewPoly()

	for i := 0; i < N; i++ {
		f.thr4FE.AggregateShares(dkShare, eknShare[i], dkShare)
	}

	f.thr4FE.AggregateShares(dkShare, eks, dkShare)

	return *dkShare

}

func (f *skEncryptionCKKS) GenerateDKShare(client4Encryption []int, yInt map[int]uint64, tsk []*drlwe.ShamirSecretShare) drlwe.ShamirSecretShare {
	yIn := make(map[int]ring.RNSScalar, len(yInt))
	for i, yi := range yInt {
		yIn[i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
	}
	dkShare := f.thr4FE.AllocateThresholdSecretShare()

	dkShare.Poly = f.ckksParams.RingQP().NewPoly()

	resultTemp := f.thr4FE.AllocateThresholdSecretShare()
	for _, pIdx := range client4Encryption {

		f.ckksParams.RingQ().MulRNSScalarMontgomery(tsk[pIdx].Q, ScalarTransform(f.ckksParams.RingQP(), yIn[pIdx]), resultTemp.Q)
		f.thr4FE.AggregateShares(dkShare, resultTemp, dkShare)
	}

	return *dkShare

}

func (f *skEncryptionCKKS) GenerateDKShareFile(client4Encryption []int, yInt map[int]uint64, id int) (drlwe.ShamirSecretShare, time.Duration) {

	start := time.Now()
	yIn := make(map[int]ring.RNSScalar, len(yInt))
	for i, yi := range yInt {
		yIn[i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
	}

	dkShare := f.thr4FE.AllocateThresholdSecretShare()

	dkShare.Poly = f.ckksParams.RingQP().NewPoly()

	resultTemp := f.thr4FE.AllocateThresholdSecretShare()

	tempEks := f.thr4FE.AllocateThresholdSecretShare()
	tempEkn := f.thr4FE.AllocateThresholdSecretShare()

	for _, pIdx := range client4Encryption {
		tempEkn.UnmarshalBinary(readSK4FEShare(id, pIdx, "n"))
		tempEks.UnmarshalBinary(readSK4FEShare(id, pIdx, "s"))
		f.ckksParams.RingQ().MulRNSScalarMontgomery(tempEks.Q, ScalarTransform(f.ckksParams.RingQP(), yIn[pIdx]), resultTemp.Q)
		f.thr4FE.AggregateShares(dkShare, resultTemp, dkShare)
		f.thr4FE.AggregateShares(dkShare, tempEkn, dkShare)
	}

	timeTotal := time.Since(start)

	start4ReadFiles := time.Now()
	readSK4FEShare(id, client4Encryption[0], "n")
	time4ReadFiles := time.Since(start4ReadFiles)
	timeClean := timeTotal - time4ReadFiles*time.Duration(2*len(client4Encryption))
	return *dkShare, timeClean

}

func (f *skEncryptionCKKS) GenerateDKShareFile4TestTime(client4Encryption []int, yInt map[int]uint64, id int) (drlwe.ShamirSecretShare, time.Duration) {

	start := time.Now()
	yIn := make(map[int]ring.RNSScalar, len(yInt))
	for i, yi := range yInt {
		yIn[i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
	}

	dkShare := f.thr4FE.AllocateThresholdSecretShare()

	dkShare.Poly = f.ckksParams.RingQP().NewPoly()

	resultTemp := f.thr4FE.AllocateThresholdSecretShare()

	tempEks := f.thr4FE.AllocateThresholdSecretShare()
	tempEkn := f.thr4FE.AllocateThresholdSecretShare()

	var time4ReadFiles time.Duration
	for i, pIdx := range client4Encryption {
		if i == 0 {
			start4ReadFiles := time.Now()
			tempEkn.UnmarshalBinary(readSK4FEShare(id, pIdx, "n"))
			tempEks.UnmarshalBinary(readSK4FEShare(id, pIdx, "s"))
			time4ReadFiles = time.Since(start4ReadFiles)
		}
		f.ckksParams.RingQ().MulRNSScalarMontgomery(tempEks.Q, ScalarTransform(f.ckksParams.RingQP(), yIn[pIdx]), resultTemp.Q)
		f.thr4FE.AggregateShares(dkShare, resultTemp, dkShare)
		f.thr4FE.AggregateShares(dkShare, tempEkn, dkShare)
	}

	timeTotal := time.Since(start)
	timeClean := timeTotal - time4ReadFiles
	return *dkShare, timeClean

}

// client4Decryption should be one-to-one correspondance to shares
// the i-th element of shares is the share of client4Decryption[i]
func (f *skEncryptionCKKS) GenerateDK(client4Decryption []int, sharesInput map[int]drlwe.ShamirSecretShare) *rlwe.SecretKey {

	points := make([]drlwe.ShamirPublicPoint, len(client4Decryption))
	for idx, pIdx := range client4Decryption {
		points[idx] = drlwe.ShamirPublicPoint(pIdx + 1)
	}

	threshold := len(points)

	cmb := NewMyCombiner(&f.ckksParams.Parameters, points, threshold)

	shares := make(map[drlwe.ShamirPublicPoint]drlwe.ShamirSecretShare)
	for i, share := range sharesInput {
		pt := drlwe.ShamirPublicPoint(i + 1)
		shares[pt] = share
	}
	result := cmb.Recover(shares)
	dk := rlwe.NewSecretKey(f.ckksParams.Parameters)
	dk.Value = result.CopyNew()
	return dk
	//fmt.Println(result.Equals(sk.Value))
}

// func (f *skEncryptionCKKS) FEDecFinalNew(dk *rlwe.SecretKey, cipherNoiseAll map[int]*rlwe.Ciphertext, cipherShareAll map[int]*rlwe.Ciphertext, yInt map[int]uint64, client4Encryption []int, result *ring.Poly) {

// 	if len(cipherShareAll) != len(yInt) {
// 		panic("The length of cipherShareAll and yInt are not match!")
// 	}

// 	y := make(map[int]ring.RNSScalar, len(yInt)+len(cipherNoiseAll))
// 	cipherAll := make(map[int]*rlwe.Ciphertext, len(yInt)+len(cipherNoiseAll))
// 	client4EncryptionExpand := make([]int, len(yInt)+len(cipherNoiseAll))

// 	count := int(0)
// 	for i, yi := range yInt {
// 		y[i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
// 		cipherAll[i] = cipherShareAll[i]
// 		client4EncryptionExpand[count] = i
// 		count = count + 1
// 	}

// 	for i, _ := range cipherNoiseAll {
// 		y[-1*i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(1)
// 		cipherAll[-1*i] = cipherNoiseAll[i]
// 		client4EncryptionExpand[count] = -1 * i
// 		count = count + 1
// 	}

// 	mesVector := f.FEDec(dk, cipherAll, y, client4EncryptionExpand)

// 	if len(result.Coeffs[0]) != len(mesVector) {
// 		panic("The length of result noise and decrypted vector are not match!")
// 	}

// 	result.Zero()

// 	for i, mi := range mesVector {
// 		result.Buff[i] = mi
// 		result.Coeffs[0][i] = mi
// 	}

//		return
//	}

func (f *skEncryptionCKKS) FEDecFinalNew(dk *rlwe.SecretKey, cipherNoiseAll map[int]*rlwe.Ciphertext, cipherShareAll map[int]*rlwe.Ciphertext, yInt map[int]uint64, client4Encryption []int, result *ring.Poly) {

	if len(cipherShareAll) != len(yInt) {
		panic("The length of cipherShareAll and yInt are not match!")
	}

	y := make(map[int]ring.RNSScalar, len(yInt))
	//cipherAll := make(map[int]*rlwe.Ciphertext, len(yInt))
	//client4EncryptionExpand := make([]int, len(yInt))

	count := int(0)
	for i, yi := range yInt {
		y[i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
		//cipherAll[i] = cipherShareAll[i]
		//client4EncryptionExpand[count] = i
		count = count + 1
	}

	mesVector := f.FEDec(dk, cipherShareAll, y, client4Encryption)

	result.Zero()
	copyLength := len(mesVector)
	if copyLength > len(result.Coeffs[0]) {
		copyLength = len(result.Coeffs[0])
	}
	for i := 0; i < copyLength; i++ {
		result.Buff[i] = mesVector[i]
		result.Coeffs[0][i] = mesVector[i]
	}

	return
}

func (f *skEncryptionCKKS) FEDecFinal(dk *rlwe.SecretKey, cipherNoiseAll map[int]*rlwe.Ciphertext, cipherShareAll map[int]*rlwe.Ciphertext, yInt map[int]uint64, client4Encryption []int, result *ring.Poly) {

	y := make(map[int]ring.RNSScalar, len(yInt)*2)
	for i, yi := range yInt {
		y[2*i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
		y[2*i+1] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(1)
	}

	cipherAll := make(map[int]*rlwe.Ciphertext, len(cipherNoiseAll)*2)
	for i, _ := range cipherNoiseAll {
		cipherAll[2*i] = cipherShareAll[i]
		cipherAll[2*i+1] = cipherNoiseAll[i]
	}

	client4EncryptionDoubleExpand := make([]int, len(client4Encryption)*2)
	for i, ci := range client4Encryption {
		client4EncryptionDoubleExpand[2*i] = ci * 2
		client4EncryptionDoubleExpand[2*i+1] = ci*2 + 1
	}

	mesVector := f.FEDec(dk, cipherAll, y, client4EncryptionDoubleExpand)

	// CKKS解码可能返回的值少于多项式系数数量
	// 只复制实际解码的值，其余设为0
	result.Zero()

	copyLength := len(mesVector)
	if copyLength > len(result.Coeffs[0]) {
		copyLength = len(result.Coeffs[0])
	}

	for i := 0; i < copyLength; i++ {
		result.Buff[i] = mesVector[i]
		result.Coeffs[0][i] = mesVector[i]
	}

	return
}

// FEDecFinalWithCombinedKey 使用组合后的密钥进行最终解密
func (f *skEncryptionCKKS) FEDecFinalWithCombinedKey(combinedKey *ring.Poly, cipherNoiseAll map[int]*rlwe.Ciphertext, cipherShareAll map[int]*rlwe.Ciphertext, yInt map[int]uint64, client4Encryption []int, result *ring.Poly) {

	y := make(map[int]ring.RNSScalar, len(yInt)*2)
	for i, yi := range yInt {
		y[2*i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
		y[2*i+1] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(1)
	}

	cipherAll := make(map[int]*rlwe.Ciphertext, len(cipherNoiseAll)*2)
	for i, _ := range cipherNoiseAll {
		cipherAll[2*i] = cipherShareAll[i]
		cipherAll[2*i+1] = cipherNoiseAll[i]
	}

	client4EncryptionDoubleExpand := make([]int, len(client4Encryption)*2)
	for i, ci := range client4Encryption {
		client4EncryptionDoubleExpand[2*i] = ci * 2
		client4EncryptionDoubleExpand[2*i+1] = ci*2 + 1
	}

	// 使用组合后的密钥创建秘密密钥对象
	combinedSK := rlwe.NewSecretKey(f.ckksParams.Parameters)
	combinedSK.Value.Q = combinedKey.CopyNew()

	mesVector := f.FEDec(combinedSK, cipherAll, y, client4EncryptionDoubleExpand)

	result.Zero()
	copyLength := len(mesVector)
	if copyLength > len(result.Coeffs[0]) {
		copyLength = len(result.Coeffs[0])
	}
	for i := 0; i < copyLength; i++ {
		result.Buff[i] = mesVector[i]
		result.Coeffs[0][i] = mesVector[i]
	}

	return
}

// FEDecFinalWithAdditiveHEKey 专门为FastShamirApproxSS处理来自加法全同态加密的组合密钥
func (f *skEncryptionCKKS) FEDecFinalWithAdditiveHEKey(combinedKey *ring.Poly, additiveHEParams ckks.Parameters, cipherNoiseAll map[int]*rlwe.Ciphertext, cipherShareAll map[int]*rlwe.Ciphertext, yInt map[int]uint64, client4Encryption []int, result *ring.Poly) {

	// 对于FastShamirApproxSS，我们需要直接解密CTsi和CTni，而不需要复杂的FE操作
	// 因为组合密钥来自加法全同态加密，参数结构不同

	y := make(map[int]ring.RNSScalar, len(yInt)*2)
	for i, yi := range yInt {
		y[2*i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
		y[2*i+1] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(1)
	}

	cipherAll := make(map[int]*rlwe.Ciphertext, len(cipherNoiseAll)*2)
	for i, _ := range cipherNoiseAll {
		cipherAll[2*i] = cipherShareAll[i]
		cipherAll[2*i+1] = cipherNoiseAll[i]
	}

	client4EncryptionDoubleExpand := make([]int, len(client4Encryption)*2)
	for i, ci := range client4Encryption {
		client4EncryptionDoubleExpand[2*i] = ci * 2
		client4EncryptionDoubleExpand[2*i+1] = ci*2 + 1
	}

	// 创建一个使用双重加密参数的临时密钥
	// 将加法全同态密钥扩展到双重加密参数的格式
	tempSK := rlwe.NewSecretKey(f.ckksParams.Parameters)

	// 复制combinedKey到双重加密参数的所有Q层级
	for i := 0; i < len(f.ckksParams.Q()); i++ {
		if i < tempSK.Value.Q.Level()+1 && i < len(combinedKey.Coeffs) {
			copy(tempSK.Value.Q.Coeffs[i], combinedKey.Coeffs[0])
		}
	}

	mesVector := f.FEDec(tempSK, cipherAll, y, client4EncryptionDoubleExpand)

	result.Zero()
	copyLength := len(mesVector)
	if copyLength > len(result.Coeffs[0]) {
		copyLength = len(result.Coeffs[0])
	}
	for i := 0; i < copyLength; i++ {
		result.Buff[i] = mesVector[i]
		result.Coeffs[0][i] = mesVector[i]
	}

	return
}

func (f *skEncryptionCKKS) FEDec(dk *rlwe.SecretKey, cipherAll map[int]*rlwe.Ciphertext, yIn map[int]ring.RNSScalar, client4Encryption []int) []uint64 {

	ringQ := f.ckksParams.RingQ()

	ct := ckks.NewCiphertext(f.ckksParams, 1, f.ckksParams.MaxLevel())

	ct.Value[1] = cipherAll[client4Encryption[0]].Value[1].CopyNew()
	ct.Value[0] = cipherAll[client4Encryption[0]].Value[0].CopyNew()

	ringQ.MulRNSScalarMontgomery(ct.Value[0], ScalarTransform(f.ckksParams.RingQP(), yIn[client4Encryption[0]]), ct.Value[0])

	ciphertextTemp := ringQ.NewPoly()

	for i, ci := range client4Encryption {
		if i != 0 {
			if cipherAll[ci].Value[1].Equals(ct.Value[1]) {
				ringQ.MulRNSScalarMontgomery(cipherAll[ci].Value[0], ScalarTransform(f.ckksParams.RingQP(), yIn[ci]), ciphertextTemp)
				ringQ.Add(ct.Value[0], ciphertextTemp, ct.Value[0])
			} else {
				fmt.Println("Wrong ciphertexts for FE: dismatch of tag polynomial for ciphertext", ci)
			}
		}
	}

	decryptor := ckks.NewDecryptor(f.ckksParams, dk)
	pt := decryptor.DecryptNew(ct)

	// CKKS解码返回complex128，需要转换为uint64
	resultComplex := f.encoder.Decode(pt, f.ckksParams.LogSlots())
	result := make([]uint64, len(resultComplex))
	for i := 0; i < len(resultComplex); i++ {
		result[i] = uint64(real(resultComplex[i]))
	}

	return result

}

// GenerateEKOnly 只生成密钥，不进行秘密分享（用于FastShamirApproxSS）
func (f *skEncryptionCKKS) GenerateEKOnly() (*rlwe.SecretKey, time.Duration) {
	start := time.Now()
	keyGen4FE := ckks.NewKeyGenerator(f.ckksParams)
	sk := keyGen4FE.GenSecretKey()
	timeTrue := time.Since(start)
	return sk, timeTrue
}

// GenerateAdditiveHEParametersPublic 生成适合加密CKKS私钥的加法全同态加密参数（公开版本用于测试）
func GenerateAdditiveHEParametersPublic(originalParams ckks.Parameters) ckks.Parameters {
	return generateAdditiveHEParameters(originalParams)
}

// generateAdditiveHEParameters 生成适合加密CKKS私钥的加法全同态加密参数
func generateAdditiveHEParameters(originalParams ckks.Parameters) ckks.Parameters {
	logN := originalParams.LogN()

	// 根据LogN选择合适的密文模数比特大小
	var minQBits int
	switch logN {
	case 10:
		minQBits = 27
	case 11:
		minQBits = 30
	case 12:
		minQBits = 33
	default:
		minQBits = 30
	}

	// 尝试不同大小的密文模数
	for qBits := minQBits; qBits <= minQBits+6; qBits++ {
		candidateQ := findNTTFriendlyPrime(logN, qBits)
		if candidateQ == 0 {
			continue
		}

		newQ := []uint64{candidateQ}

		// 构建RLWE参数
		heParams, err := rlwe.NewParameters(
			originalParams.LogN(),
			newQ,
			nil,
			0,
			originalParams.HammingWeight(),
			originalParams.Sigma(),
			originalParams.RingType(),
			rlwe.NewScale(1),
			originalParams.DefaultNTTFlag(),
		)

		if err != nil {
			continue
		}

		// 构建CKKS参数
		ckksParams, err := ckks.NewParameters(heParams, originalParams.LogSlots())
		if err != nil {
			continue
		}

		return ckksParams
	}

	panic("Cannot find suitable NTT-friendly parameters for additive homomorphic encryption!")
}

// generatePrimeOfBitSize 生成指定比特大小的素数
func generatePrimeOfBitSize(bitSize int) uint64 {
	if bitSize <= 0 || bitSize > 63 {
		fmt.Printf("[DEBUG] Invalid bitSize: %d\n", bitSize)
		return 0
	}

	// 生成指定比特大小范围内的候选数
	min := uint64(1) << (bitSize - 1)
	max := (uint64(1) << bitSize) - 1

	fmt.Printf("[DEBUG] generatePrimeOfBitSize: bitSize=%d, min=%d, max=%d\n", bitSize, min, max)

	// 限制搜索范围以避免无限循环
	maxAttempts := 100 // 减少尝试次数

	for i := 0; i < maxAttempts; i++ {
		// 更简单的候选数生成
		candidate := min + uint64(i*997) // 使用质数步长
		if candidate > max {
			candidate = min + uint64(i*2+1) // 确保奇数
		}

		// 确保是奇数
		if candidate%2 == 0 {
			candidate++
		}

		// 进行素性测试
		if big.NewInt(int64(candidate)).ProbablyPrime(10) {
			fmt.Printf("[DEBUG] Found prime: %d\n", candidate)
			return candidate
		}
	}

	fmt.Printf("[DEBUG] Failed to find prime of bitSize %d\n", bitSize)
	return 0
}

// findNTTFriendlyPrime 找到满足 q ≡ 1 (mod 2N) 且为素数的模数
func findNTTFriendlyPrime(logN int, minBits int) uint64 {
	N := uint64(1 << logN)
	modulus := 2 * N

	// 从合适的起始点开始搜索
	start := uint64(1) << minBits
	// 确保起始点满足模数条件
	if start%modulus != 1 {
		start = start + (modulus - (start % modulus)) + 1
	}

	for candidate := start; candidate < start+modulus*1000; candidate += modulus {
		if big.NewInt(int64(candidate)).ProbablyPrime(10) {
			return candidate
		}
	}
	return 0
}

// FEDecFinalWithECElGamalKey 专门为Fast2ShamirApproxSS处理来自EC-ElGamal门限解密的组合密钥
func (f *skEncryptionCKKS) FEDecFinalWithECElGamalKey(combinedKey *ring.Poly, originalParams ckks.Parameters, cipherNoiseAll map[int]*rlwe.Ciphertext, cipherShareAll map[int]*rlwe.Ciphertext, yInt map[int]uint64, client4Encryption []int, result *ring.Poly) {

	// 对于Fast2ShamirApproxSS，我们需要直接解密CTsi和CTni
	// 组合密钥来自EC-ElGamal门限解密

	y := make(map[int]ring.RNSScalar, len(yInt)*2)
	for i, yi := range yInt {
		y[2*i] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(yi)
		y[2*i+1] = f.ckksParams.RingQP().NewRNSScalarFromUInt64(1)
	}

	cipherAll := make(map[int]*rlwe.Ciphertext, len(cipherNoiseAll)*2)
	for i, _ := range cipherNoiseAll {
		cipherAll[2*i] = cipherShareAll[i]
		cipherAll[2*i+1] = cipherNoiseAll[i]
	}

	client4EncryptionDoubleExpand := make([]int, len(client4Encryption)*2)
	for i, ci := range client4Encryption {
		client4EncryptionDoubleExpand[2*i] = ci * 2
		client4EncryptionDoubleExpand[2*i+1] = ci*2 + 1
	}

	// 创建一个使用双重加密参数的临时密钥
	// 将EC-ElGamal解密得到的密钥扩展到双重加密参数的格式
	tempSK := rlwe.NewSecretKey(f.ckksParams.Parameters)

	// 复制combinedKey到双重加密参数的所有Q层级
	for i := 0; i < len(f.ckksParams.Q()); i++ {
		if i < tempSK.Value.Q.Level()+1 && i < len(combinedKey.Coeffs) {
			copy(tempSK.Value.Q.Coeffs[i], combinedKey.Coeffs[0])
		}
	}

	mesVector := f.FEDec(tempSK, cipherAll, y, client4EncryptionDoubleExpand)

	result.Zero()
	copyLength := len(mesVector)
	if copyLength > len(result.Coeffs[0]) {
		copyLength = len(result.Coeffs[0])
	}
	for i := 0; i < copyLength; i++ {
		result.Buff[i] = mesVector[i]
		result.Coeffs[0][i] = mesVector[i]
	}

	return
}
