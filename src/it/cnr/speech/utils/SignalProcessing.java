package it.cnr.speech.utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import it.cnr.speech.filters.LowPassFilterDynamic;

/**
 * includes tools for basic signal transformations: delta + double delta center frequency cepstral coefficients calculation spectrum frequency cut transformation to and from Rapid Miner Example Set filterbanks fequency to mel frequency to index in fft sinusoid signal generation inverse mel log10 mel filterbanks sample to time and time to sample signal timeline generation index to time in spectrogram spectrogram calculation and display time to index in spectrogram
 * 
 * @author coro
 * 
 */
public class SignalProcessing {

	public int windowShiftSamples;
	public int windowSizeSamples;
	public int samplingRate;
	public double signal[];
	
	public static int timeToSamples(double time, double fs) {
		return (int) Math.round(fs * time);
	}
	public static int frequencyIndex(float frequency, int fftSize, float samplingRate) {
		return Math.round(frequency * fftSize / samplingRate);
	}
	public void getSignal(File audio) throws Exception{
		AudioBits bits = new AudioBits(audio);
		 signal = bits.getDoubleVectorAudio();
		samplingRate = (int) bits.getAudioFormat().getSampleRate(); 
		bits.ais.close();
	}
	
	/**
	 * Downsamples the loaded signal using windowed-sinc interpolation
	 * with an anti-aliasing low-pass filter.
	 *
	 * Supports integer and non-integer sampling-rate ratios.
	 */
	public void downsampleSignal(int targetSamplingRate) {
	    if (signal == null || signal.length == 0) {
	        throw new IllegalStateException("Load a non-empty signal first.");
	    }
	    if (samplingRate <= 0 || targetSamplingRate <= 0) {
	        throw new IllegalArgumentException("Sampling rates must be positive.");
	    }
	    if (targetSamplingRate > samplingRate) {
	        throw new IllegalArgumentException(
	            "Target sampling rate must not exceed the original rate."
	        );
	    }
	    if (targetSamplingRate == samplingRate) {
	        return;
	    }

	    final double ratio = (double) targetSamplingRate / samplingRate;

	    // Cutoff in cycles per input sample:
	    // 90% of the target Nyquist frequency leaves a transition band.
	    final double cutoff = 0.45 * ratio;

	    // Increase filter length for larger downsampling factors.
	    final int radius = (int) Math.ceil(32.0 / ratio);

	    final int outputLength = (int) Math.ceil(signal.length * ratio);
	    final double[] resampled = new double[outputLength];

	    for (int i = 0; i < outputLength; i++) {
	        final double position = i / ratio;

	        final int start = (int) Math.max(
	            0.0, Math.ceil(position - radius)
	        );
	        final int end = (int) Math.min(
	            signal.length - 1.0, Math.floor(position + radius)
	        );

	        double weightedSum = 0.0;
	        double weightSum = 0.0;

	        for (int j = start; j <= end; j++) {
	            final double distance = position - j;
	            final double x = 2.0 * cutoff * distance;

	            final double sinc = Math.abs(x) < 1e-12
	                ? 1.0
	                : Math.sin(Math.PI * x) / (Math.PI * x);

	            // Symmetric Hann window.
	            final double window = 0.5 * (
	                1.0 + Math.cos(Math.PI * distance / radius)
	            );

	            final double weight = 2.0 * cutoff * sinc * window;

	            weightedSum += signal[j] * weight;
	            weightSum += weight;
	        }

	        // Preserve constant-signal amplitude, including at the boundaries.
	        resampled[i] = weightedSum / weightSum;
	    }

	    signal = resampled;
	    samplingRate = targetSamplingRate;
	}
	
	public static double calculateAverageEnvelopeLevel(double[] signal) {
        double sum = 0.0;

        for (double sample : signal) {
            sum += Math.abs(sample);
        }

        return sum / signal.length;
    }
	
	public static double[] normaliseByAverageEnvelopeLevel (double signal[]) {
		
		double ael = calculateAverageEnvelopeLevel(signal);

		if (!Double.isFinite(ael)) {
		    throw new IllegalArgumentException("Non-finite envelope level.");
		}

		double[] result = new double[signal.length];

		if (ael == 0.0) {
		    return result;
		}

		for (int i = 0; i < signal.length; i++) {
		    result[i] = signal[i] / ael;
		}
		return result;
		/*
		double ael = calculateAverageEnvelopeLevel(signal);
		double aelSignal [] = new double[signal.length];
		for (int i=0;i<signal.length;i++) {
			aelSignal[i] = signal[i]/ael;
		}
		return aelSignal;
		*/
	}
	
	public static double samplesToTime(int samples, double fs) {
		return (double) samples / fs;
	}
	
	/*
	public double[][] shortTermFFT(File audio, double windowSize, double windowShift) throws Exception{
		
		getSignal(audio);
		return shortTermFFT(signal, samplingRate ,windowSize, windowShift);
	}
	*/
	public double[][] shortTermFFT(double[] signal, int samplingRate ,double windowSize, double windowShift) throws Exception{	
		windowSizeSamples = SignalProcessing.timeToSamples(windowSize, samplingRate);
		System.out.println("Original window sample: "+windowSize+"s"+" "+windowSizeSamples+" (samples)");
		windowSizeSamples = UtilsMath.powerTwoApproximation(windowSizeSamples);
		System.out.println("Approx window sample: "+SignalProcessing.samplesToTime(windowSizeSamples,samplingRate) +"s"+" "+windowSizeSamples+" (samples)");
		
		windowShiftSamples = SignalProcessing.timeToSamples(windowShift, samplingRate);
		windowShiftSamples = UtilsMath.powerTwoApproximation(windowShiftSamples);
		
		System.out.println("Running FFT with "+windowSizeSamples+" by "+windowShiftSamples+" ...");
		
		List<double[]> spectra = new ArrayList<double[]>();
        
        for (int i = 0; i < signal.length; i += windowShiftSamples) {
        	if ((i+windowSizeSamples)>signal.length)
        		break;
        	
            // Extract a windowed segment of the signal
            double[] windowedSegment = LowPassFilterDynamic.getWindowedSegment(signal, i, windowSizeSamples);
            windowedSegment = LowPassFilterDynamic.hammingWindow(windowedSegment);
            
            // Compute the Fourier Transform for the windowed segment
            Complex[] complexSpectrum = LowPassFilterDynamic.computeFourierTransform(windowedSegment);

            // Apply bandpass filter to the spectrum
            double[] absSpectrum = new double[complexSpectrum.length];
            for (int k=0;k<absSpectrum.length;k++) {
            	
            	absSpectrum[k] = complexSpectrum[k].abs();
            }
            
            spectra.add(absSpectrum);
        }

        double[][] spectrum = spectra.toArray(new double[spectra.size()][]);
        
        return spectrum;
 
	}
	
	public static double[][] cutSpectrum(double[][] spectrum, float minFreq, float maxfreq, int fftWindowSize, int samplingRate) {
		int minFrequencyIndex = frequencyIndex(minFreq, fftWindowSize, samplingRate);
		int maxFrequencyIndex = frequencyIndex(maxfreq, fftWindowSize, samplingRate);

		double[][] cutSpectrum = new double[spectrum.length][maxFrequencyIndex - minFrequencyIndex + 1];

		for (int i = 0; i < spectrum.length; i++) {
			cutSpectrum[i] = Arrays.copyOfRange(spectrum[i], minFrequencyIndex, maxFrequencyIndex+1);
		}

		return cutSpectrum;
	}
	
	
	/**
	 * Filters an envelope with a 16th-order Butterworth low-pass.
	 *
	 * The cutoff is the half-power frequency (approximately -3 dB).
	 * Uses eight cascaded biquads for numerical stability.
	 *
	 * Returns a new array without modifying the input.
	 */
	public static double[] lowPassEnvelope(
	        double[] input, double samplingRate, double cutoffHz) {

	    if (input == null) {
	        throw new IllegalArgumentException("Input must not be null.");
	    }
	    if (!Double.isFinite(samplingRate)
	            || !Double.isFinite(cutoffHz)
	            || samplingRate <= 0
	            || cutoffHz <= 0
	            || cutoffHz >= samplingRate / 2.0) {
	        throw new IllegalArgumentException(
	            "Cutoff must be positive and below Nyquist."
	        );
	    }

	    final int order = 16;
	    double[] output = input.clone();

	    final double omega = 2.0 * Math.PI * cutoffHz / samplingRate;
	    final double sinOmega = Math.sin(omega);
	    final double cosOmega = Math.cos(omega);

	    // Equivalent to (1 - cos(omega)) / 2, with better
	    // numerical accuracy at low cutoff frequencies.
	    final double sinHalf = Math.sin(omega / 2.0);
	    final double numerator = sinHalf * sinHalf;

	    // Lower-Q sections first.
	    for (int section = 0; section < order / 2; section++) {
	        final double q = 1.0 / (
	            2.0 * Math.cos(
	                (2.0 * section + 1.0) * Math.PI / (2.0 * order)
	            )
	        );

	        final double alpha = sinOmega / (2.0 * q);
	        final double a0 = 1.0 + alpha;

	        final double b0 = numerator / a0;
	        final double b1 = 2.0 * b0;
	        final double b2 = b0;
	        final double a1 = -2.0 * cosOmega / a0;
	        final double a2 = (1.0 - alpha) / a0;

	        // Transposed direct-form II, initially at rest.
	        double state1 = 0.0;
	        double state2 = 0.0;

	        for (int n = 0; n < output.length; n++) {
	            final double x = output[n];
	            final double y = b0 * x + state1;

	            state1 = b1 * x - a1 * y + state2;
	            state2 = b2 * x - a2 * y;

	            output[n] = y;
	        }
	    }

	    return output;
	}
	
	
	/**
	 * Calculates the magnitude of the 4 Hz DFT coefficient
	 * over 20-sample Hamming windows with a one-sample hop.
	 *
	 * Input must be sampled at 80 Hz.
	 * Only complete windows are included.
	 */
	public static double[] modulationMagnitude4Hz(
	        double[] envelope, int samplingRate) {

	    if (envelope == null) {
	        throw new IllegalArgumentException("Envelope must not be null.");
	    }
	    if (samplingRate != 80) {
	        throw new IllegalArgumentException(
	            "This method requires an envelope sampled at 80 Hz."
	        );
	    }

	    final int windowSamples = 20;
	    final int hopSamples = 1;
	    final double frequencyHz = 4.0;

	    if (envelope.length < windowSamples) {
	        throw new IllegalArgumentException(
	            "At least 20 envelope samples (250 ms) are required."
	        );
	    }

	    final int frameCount =
	        1 + (envelope.length - windowSamples) / hopSamples;

	    final double[] magnitudes = new double[frameCount];

	    // Precompute the windowed DFT basis.
	    final double[] realWeights = new double[windowSamples];
	    final double[] imagWeights = new double[windowSamples];

	    for (int n = 0; n < windowSamples; n++) {
	        final double hamming =
	            0.54 - 0.46 * Math.cos(
	                2.0 * Math.PI * n / (windowSamples - 1)
	            );

	        final double angle =
	            2.0 * Math.PI * frequencyHz * n / samplingRate;

	        realWeights[n] = hamming * Math.cos(angle);
	        imagWeights[n] = -hamming * Math.sin(angle);
	    }

	    for (int frame = 0; frame < frameCount; frame++) {
	        final int start = frame * hopSamples;

	        double real = 0.0;
	        double imaginary = 0.0;

	        for (int n = 0; n < windowSamples; n++) {
	            final double sample = envelope[start + n];

	            real += sample * realWeights[n];
	            imaginary += sample * imagWeights[n];
	        }

	        // Unnormalized magnitude, matching the existing forward FFT.
	        magnitudes[frame] = Math.hypot(real, imaginary);
	    }

	    return magnitudes;
	}
}
