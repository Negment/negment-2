// noise-processor.js
class NoiseProcessor extends AudioWorkletProcessor {
    constructor() {
        super();
        this.fftSize = 4096;
        this.fft = this.createFFT(this.fftSize);
        this.re = new Float32Array(this.fftSize);
        this.im = new Float32Array(this.fftSize).fill(0);
        this.bufferIndex = 0;
    }
    
    static get parameterDescriptors() {
        return [
            { name: 'minFrequency', defaultValue: 0 },
            { name: 'maxFrequency', defaultValue: 20000 },
            { name: 'gain', defaultValue: 0.1 }
        ];
    }
    
    // FFTクラスをWorklet内に定義
    createFFT(size) {
        function FFT(size) {
            this.size = size;
            this.size_log = Math.log2(size);
        }
    
        FFT.prototype.transform = function(re, im) {
            let m = this.size_log;
            let n = this.size;
            
            for (let k = 0; k < n; k++) {
                let j = 0;
                for (let i = 0; i < m; i++) {
                    j |= (((k >> i) & 1) << (m - 1 - i));
                }
                if (j > k) {
                    [re[k], re[j]] = [re[j], re[k]];
                    [im[k], im[j]] = [im[j], im[k]];
                }
            }
            
            let a, b, c, s, angle, len_csize, len_csize_s, tr;
            
            for (let h = 2; h <= n; h <<= 1) {
                angle = -2 * Math.PI / h;
                len_csize_s = Math.cos(angle);
                len_csize = Math.sin(angle);
                a = 1;
                b = 0;
                for (let j = 0; j < h / 2; j++) {
                    for (let k = j; k < n; k += h) {
                        c = re[k + h / 2] * a - im[k + h / 2] * b;
                        s = re[k + h / 2] * b + im[k + h / 2] * a;
                        re[k + h / 2] = re[k] - c;
                        im[k + h / 2] = im[k] - s;
                        re[k] += c;
                        im[k] += s;
                    }
                    tr = a;
                    a = a * len_csize_s - b * len_csize;
                    b = tr * len_csize + b * len_csize_s;
                }
            }
        };
        
        FFT.prototype.inverseTransform = function(re, im) {
            let n = this.size;
            im = im.map(val => -val);
            this.transform(re, im);
            im = im.map(val => -val);
            re = re.map(val => val / n);
            im = im.map(val => val / n);
            return [re, im];
        };
    
        return new FFT(size);
    }
    
    // RMSパワーを計算
    calculateRMS(buffer) {
        let sumOfSquares = 0;
        for (let i = 0; i < buffer.length; i++) {
            sumOfSquares += buffer[i] * buffer[i];
        }
        return Math.sqrt(sumOfSquares / buffer.length);
    }
    
    process(inputs, outputs, parameters) {
        const output = outputs[0][0];
        const minFrequency = parameters.minFrequency[0];
        const maxFrequency = parameters.maxFrequency[0];
        const gainValue = parameters.gain[0];
        const sampleRate = 44100; // AudioWorkletでは固定値

        // リアルタイムでノイズを生成
        for (let i = 0; i < this.re.length; i++) {
            this.re[i] = Math.random() * 2 - 1;
            this.im[i] = 0;
        }

        // FFTと周波数カット
        this.fft.transform(this.re, this.im);
        const freqResolution = sampleRate / this.fftSize;
        const minBin = Math.floor(minFrequency / freqResolution);
        const maxBin = Math.floor(maxFrequency / freqResolution);

        for (let i = 0; i < this.fftSize / 2; i++) {
            if (i < minBin || i > maxBin) {
                this.re[i] = 0;
                this.im[i] = 0;
                this.re[this.fftSize - i] = 0;
                this.im[this.fftSize - i] = 0;
            }
        }
        
        // IFFTと音量補正
        const [outRe] = this.fft.inverseTransform(this.re, this.im);
        
        let maxAmp = 0;
        for (let i = 0; i < outRe.length; i++) {
            if (Math.abs(outRe[i]) > maxAmp) {
                maxAmp = Math.abs(outRe[i]);
            }
        }
        
        if (maxAmp > 0) {
            const normalizedNoise = new Float32Array(outRe.length);
            for (let i = 0; i < outRe.length; i++) {
                normalizedNoise[i] = outRe[i] / maxAmp;
            }
            
            const rmsPower = this.calculateRMS(normalizedNoise);
            const correctedGain = (rmsPower > 0) ? gainValue / rmsPower : 0;
            
            for (let i = 0; i < output.length; i++) {
                output[i] = normalizedNoise[i] * correctedGain;
            }
        } else {
            for (let i = 0; i < output.length; i++) {
                output[i] = 0;
            }
        }
        
        return true;
    }
}

registerProcessor('noise-processor', NoiseProcessor);
