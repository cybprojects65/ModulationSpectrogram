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
	
	public static double calculateAverageEnvelopeLevel(double[] signal) {
        double sum = 0.0;

        for (double sample : signal) {
            sum += Math.abs(sample);
        }

        return sum / signal.length;
    }
	
	public static double[] normaliseByAverageEnvelopeLevel (double signal[]) {
		
		double ael = calculateAverageEnvelopeLevel(signal);
		double aelSignal [] = new double[signal.length];
		for (int i=0;i<signal.length;i++) {
			aelSignal[i] = signal[i]/ael;
		}
		return aelSignal;
	}
	
	public static double samplesToTime(int samples, double fs) {
		return (double) samples / fs;
	}
	
	public double[][] shortTermFFT(File audio, double windowSize, double windowShift) throws Exception{
		
		getSignal(audio);
		return shortTermFFT(signal, samplingRate ,windowSize, windowShift);
	}
	
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
}
