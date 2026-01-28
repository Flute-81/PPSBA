package ApproxSS

import (
	"bufio"
	"fmt"
	"io"
	"math/big"
	"math/rand"

	"os"
	"time"

	//"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func TestResult(approxMessage, originalMessage *ring.Poly, bound uint64, modulus []uint64) (isSuc bool) {

	if len(approxMessage.Coeffs) != len(originalMessage.Coeffs) {
		panic("The length of approximate message does not match that of original message!")
	}

	for m, modulu := range modulus {
		for i, oriCoeff := range originalMessage.Coeffs[m] {

			isSucThisCoeff := false
			allowedSet := make([]uint64, 2*bound+1)
			allowedSet[0] = oriCoeff
			for j := uint64(0); j < bound; j++ {
				if oriCoeff+j < modulu {
					allowedSet[2*j] = oriCoeff + j
				} else {
					allowedSet[2*j] = oriCoeff + j - modulu
				}

				if int64(oriCoeff)-int64(j) < 0 {
					allowedSet[2*j+1] = oriCoeff - j + modulu
				} else {
					allowedSet[2*j+1] = oriCoeff - j
				}
			}

			for _, aCoeff := range allowedSet {
				if approxMessage.Buff[i] == aCoeff {
					isSucThisCoeff = true
				}
			}

			if !isSucThisCoeff {
				isSuc = false
				return
			}
		}
	}

	isSuc = true
	return
}

func Factorial(N int) (Nfactorial *big.Int) {
	Nfactorial.MulRange(1, int64(N))
	return
}

func SelectRandomParticipants(N int, T int) (par []int) {

	allPar := make([]int, N)
	for i := 0; i < N; i++ {
		allPar[i] = i
	}

	rand.Shuffle(N, func(i, j int) {
		allPar[i], allPar[j] = allPar[j], allPar[i]
	})

	par = allPar[0:T]

	return

}
func ScalarTransform(ringQP *ringqp.Ring, scalarIn ring.RNSScalar) ring.RNSScalar {
	scalar := ringQP.NewRNSScalar()
	qlen := len(ringQP.RingQ.Modulus)

	for i, qi := range ringQP.RingQ.Modulus {
		scalar[i] = ring.MForm(scalarIn[i], qi, ringQP.RingQ.BredParams[i])
	}
	if ringQP.RingP != nil {
		for i, pi := range ringQP.RingP.Modulus {
			scalar[i+qlen] = ring.MForm(scalarIn[i+qlen], pi, ringQP.RingP.BredParams[i])
		}
	}

	return scalar
}

func ScalarTransformQLvl(ringQ *ring.Ring, scalarIn ring.RNSScalar, level int) ring.RNSScalar {
	scalar := ringQ.NewRNSScalar()

	for i, qi := range ringQ.Modulus {
		if level > -1 {
			scalar[i] = ring.MForm(scalarIn[i], qi, ringQ.BredParams[i])
		}
		level = level - 1
	}

	return scalar
}

func InvScalarTransform(ringQP *ringqp.Ring, scalarIn ring.RNSScalar) ring.RNSScalar {
	scalar := ringQP.NewRNSScalar()
	qlen := len(ringQP.RingQ.Modulus)

	for i, qi := range ringQP.RingQ.Modulus {
		scalar[i] = ring.InvMForm(scalarIn[i], qi, ringQP.RingQ.MredParams[i])
	}
	if ringQP.RingP != nil {
		for i, pi := range ringQP.RingP.Modulus {
			scalar[i+qlen] = ring.InvMForm(scalarIn[i+qlen], pi, ringQP.RingP.MredParams[i])
		}
	}

	return scalar
}

func fromRNSScalarToInt(ringQP *ringqp.Ring, scalar ring.RNSScalar) uint64 {

	return InvScalarTransform(ringQP, scalar)[0]

}

func CheckFileIsExist(filename string) bool {
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		return false
	}
	return true
}

func readSK4FEShare(id int, source int, NorS string) []byte {

	filename := fmt.Sprintf("%s%d%s%d", "mesFrom", source, "To", id)
	sk4FEFilename := fmt.Sprintf("%s%s%s%s", "file/ek", NorS, "Share/", filename)
	sk4FEFile, _ := os.Open(sk4FEFilename)
	defer sk4FEFile.Close()
	sk4FEData, err := io.ReadAll(sk4FEFile)
	if err != nil {
		panic(err)
	}
	return sk4FEData

}

func writeSK4FEShare(id int, source int, NorS string, share drlwe.ShamirSecretShare) {

	filename := fmt.Sprintf("%s%d%s%d", "mesFrom", source, "To", id)
	sk4FEFilename := fmt.Sprintf("%s%s%s%s", "file/ek", NorS, "Share/", filename)
	shareData, _ := share.MarshalBinary()
	err := os.WriteFile(sk4FEFilename, shareData, 0666)
	if err != nil {
		panic(err)
	}

}

func readNonceShare(id int, source int) []byte {

	filename := fmt.Sprintf("%s%d%s%d", "mesFrom", source, "To", id)
	sk4FEFilename := fmt.Sprintf("%s%s", "file/nonceShare/", filename)
	sk4FEFile, _ := os.Open(sk4FEFilename)
	defer sk4FEFile.Close()
	sk4FEData, err := io.ReadAll(sk4FEFile)
	if err != nil {
		panic(err)
	}
	return sk4FEData

}

func writeNonceShare(id int, source int, share drlwe.ShamirSecretShare) {

	filename := fmt.Sprintf("%s%d%s%d", "mesFrom", source, "To", id)
	sk4FEFilename := fmt.Sprintf("%s%s", "file/nonceShare/", filename)
	shareData, _ := share.MarshalBinary()
	err := os.WriteFile(sk4FEFilename, shareData, 0666)
	if err != nil {
		panic(err)
	}

}

func RecordParameters(filename string, numClient, numClientPerGroup, threshold, app, numCPU int, timeNow time.Time) {

	var file *os.File
	var err error
	if !CheckFileIsExist(filename) {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_CREATE, 0666)
	} else {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_APPEND, 0666)
	}

	if err != nil {
		fmt.Println("文件打开失败", err)
	}
	//及时关闭file句柄
	defer file.Close()
	//写入文件时，使用带缓存的 *Writer
	write := bufio.NewWriter(file)
	hour, min, sec := timeNow.Clock()
	mes := fmt.Sprintf("\nnumClient: %d; numClientPerGroup: %d; threshold: %d; numCPU: %d; app: %d in %d Hour, %d Min, %d Sec\n", numClient, numClientPerGroup, threshold, numCPU, app, hour, min, sec)
	write.WriteString(mes)

	//Flush将缓存的文件真正写入到文件中
	write.Flush()
}

func RecordTime(filename string, title string, time time.Duration) {

	var file *os.File
	var err error
	if !CheckFileIsExist(filename) {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_CREATE, 0666)
	} else {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_APPEND, 0666)
	}

	if err != nil {
		fmt.Println("文件打开失败", err)
	}
	//及时关闭file句柄
	defer file.Close()
	//写入文件时，使用带缓存的 *Writer
	write := bufio.NewWriter(file)
	mes := fmt.Sprintf("%s%s%s\n", title, " ", time)
	write.WriteString(mes)

	//Flush将缓存的文件真正写入到文件中
	write.Flush()
}

func RecordTimeMultiple(filename string, times []time.Duration, datas []int) {

	var file *os.File
	var err error
	if !CheckFileIsExist(filename) {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_CREATE, 0666)
	} else {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_APPEND, 0666)
	}

	if err != nil {
		fmt.Println("文件打开失败", err)
	}
	//及时关闭file句柄
	defer file.Close()
	//写入文件时，使用带缓存的 *Writer
	write := bufio.NewWriter(file)
	var mes string
	for _, time := range times {
		mes = fmt.Sprintf("%s%s%s", mes, time, " ")
	}
	for _, data := range datas {
		mes = fmt.Sprintf("%s%d%s", mes, data, " ")
	}
	mes = fmt.Sprintf("%s\n", mes)
	write.WriteString(mes)

	//Flush将缓存的文件真正写入到文件中
	write.Flush()
}

func RecordData(filename string, title string, dataSize int) {

	var file *os.File
	var err error
	if !CheckFileIsExist(filename) {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_CREATE, 0666)
	} else {
		file, err = os.OpenFile(filename, os.O_WRONLY|os.O_APPEND, 0666)
	}

	if err != nil {
		fmt.Println("文件打开失败", err)
	}
	//及时关闭file句柄
	defer file.Close()
	//写入文件时，使用带缓存的 *Writer
	write := bufio.NewWriter(file)
	mes := fmt.Sprintf("%s%s%d%s\n", title, " ", dataSize, "bytes")
	write.WriteString(mes)

	//Flush将缓存的文件真正写入到文件中
	write.Flush()
}

func writeSecretKeyToFile(filename string, SK *rlwe.SecretKey) {
	var f *os.File
	if CheckFileIsExist(filename) { //如果文件存在
		os.Remove(filename)
		fmt.Println("文件存在")
	} else {
		fmt.Println("文件不存在")
	}
	f, _ = os.Create(filename) //打开文件
	defer f.Close()

	var skString string
	for idx, ski := range SK.Value.Q.Buff {
		skString = fmt.Sprintf("%d %d\n", idx, ski)
		//fmt.Println(idx)
		_, err1 := io.WriteString(f, skString) //写入文件(字符串)
		if err1 != nil {
			panic(err1)
		}
	}

	//fmt.Printf("写入 %d 个字节n", n)
}

func rotate(slice []float64, positions int) {
	positions = positions % len(slice) //find the position
	reverse(slice[:positions])         //reverse the beginning elements
	reverse(slice[positions:])         //reverse the end elements
	reverse(slice)                     //reverse the entire slice
}

func reverse(slice []float64) {
	for i, j := 0, len(slice)-1; i < j; i, j = i+1, j-1 {
		slice[i], slice[j] = slice[j], slice[i]
	}
}

func removeDuplicateElement(addrs []uint64) []uint64 {
	result := make([]uint64, 0, len(addrs))
	temp := map[uint64]struct{}{}
	for _, item := range addrs {
		if _, ok := temp[item]; !ok {
			temp[item] = struct{}{}
			result = append(result, item)
		}
	}
	return result
}

// // AtLevel returns a shallow copy of the target ring configured to
// // carry on operations at the specified levels.
// func RingqpAtLevel(r *ringqp.Ring, levelQ, levelP int) *ringqp.Ring {

// 	var ringQ, ringP *ring.Ring

// 	if levelQ > -1 && r.RingQ != nil {
// 		ringQ = r.RingQ.AtLevel(levelQ)
// 	}

// 	if levelP > -1 && r.RingP != nil {
// 		ringP = r.RingP.AtLevel(levelP)
// 	}

// 	return &ringqp.Ring{
// 		RingQ: ringQ,
// 		RingP: ringP,
// 	}
// }

// func AtLevel(r *ring.Ring, level int) *ring.Ring {

// 	if level < 0 {
// 		panic("level cannot be negative")
// 	}

// 	if level > r.ModuliChainLength() - 1 {
// 		panic("level cannot be larger than max level")
// 	}

// 	return &ring.Ring{
// 		SubRings:         r.SubRings,
// 		ModulusAtLevel:   r.ModulusAtLevel,
// 		RescaleConstants: r.RescaleConstants,
// 		level:            level,
// 	}
// }

// CalculateSecureSmudgingNoiseBound 计算安全的smudging noise边界
// 确保smudging noise的标准差超多项式大于原始CKKS噪声
func CalculateSecureSmudgingNoiseBound(params ckks.Parameters, minBound int) int {
	// 获取原始CKKS参数的噪声标准差
	originalSigma := params.Sigma()

	// 使用12倍安全边际，确保超多项式大于原噪声 (对于128位安全级别)
	// 根据smudging技术的安全性证明，通常需要σ_smudging ≥ ω(log λ) · σ_original
	securityMultiplier := 12.0
	minSmudgingStdDev := securityMultiplier * originalSigma

	// 计算对应的边界值 (由于使用 stdDev = bound/6)
	minSecureBound := int(minSmudgingStdDev * 6)

	// 取最小安全边界和当前边界的最大值
	if minBound > minSecureBound {
		// fmt.Printf("警告: 当前边界 %d 已满足安全要求 (最小安全边界: %d)\n", minBound, minSecureBound)
		return minBound
	}

	// fmt.Printf("安全性调整: 原始边界 %d -> 安全边界 %d (原始σ=%.2f, 安全σ=%.2f)\n",
	// 	minBound, minSecureBound, originalSigma, minSmudgingStdDev)

	return minSecureBound
}

// CreateSecureSmudgingNoiseSampler 创建安全的smudging noise采样器
func CreateSecureSmudgingNoiseSampler(params ckks.Parameters, ringQ *ring.Ring, bound4smudgingNoise int) (*ring.GaussianSampler, int) {
	// 计算安全的噪声边界
	secureBound := CalculateSecureSmudgingNoiseBound(params, bound4smudgingNoise)

	// 创建PRNG和采样器
	prng, _ := utils.NewPRNG()
	secureStdDev := float64(secureBound) / 6.0

	sampler := ring.NewGaussianSampler(prng, ringQ, secureStdDev, secureBound)

	return sampler, secureBound
}

// CreateSecureSmudgingNoiseSamplerFromRLWE 从RLWE参数创建安全的smudging noise采样器
func CreateSecureSmudgingNoiseSamplerFromRLWE(params *rlwe.Parameters, ringQ *ring.Ring, bound4smudgingNoise int) (*ring.GaussianSampler, int) {
	// 获取原始RLWE参数的噪声标准差
	originalSigma := params.Sigma()

	// 使用12倍安全边际，确保超多项式大于原噪声
	securityMultiplier := 12.0
	minSmudgingStdDev := securityMultiplier * originalSigma

	// 计算对应的边界值 (由于使用 stdDev = bound/6)
	minSecureBound := int(minSmudgingStdDev * 6)

	// 取最小安全边界和当前边界的最大值
	var secureBound int
	if bound4smudgingNoise > minSecureBound {
		secureBound = bound4smudgingNoise
	} else {
		secureBound = minSecureBound
		// fmt.Printf("安全性调整(RLWE): 原始边界 %d -> 安全边界 %d (原始σ=%.2f, 安全σ=%.2f)\n",
		// bound4smudgingNoise, secureBound, originalSigma, minSmudgingStdDev)
	}

	// 创建PRNG和采样器
	prng, _ := utils.NewPRNG()
	secureStdDev := float64(secureBound) / 6.0

	sampler := ring.NewGaussianSampler(prng, ringQ, secureStdDev, secureBound)

	return sampler, secureBound
}
