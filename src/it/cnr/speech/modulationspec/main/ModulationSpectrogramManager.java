package it.cnr.speech.modulationspec.main;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.speech.utils.UtilsMath;

public class ModulationSpectrogramManager {

	ModulationSpectrogram ms = new ModulationSpectrogram();

	public List<double[]> extractModulationSpectrogram(File entireSignal, int modulation_spectrogram_nfeatures,
			double modulation_spectrogram_max_frequency, boolean delta, boolean doubledelta) throws Exception {

		List<double[]> features = new ArrayList<double[]>();

		boolean saturate = true;
		File output = null;
		boolean addDeltas = false;

		System.out.println("Calculating Modulation Spectrogram with " + modulation_spectrogram_nfeatures
				+ " features and a max frequency of " + modulation_spectrogram_max_frequency);
		ms.calcMS(entireSignal, output, saturate, addDeltas, modulation_spectrogram_nfeatures,
				modulation_spectrogram_max_frequency);
		double[][] mspec = UtilsMath.traspose(ms.modulationSpectrogram);
		int nfeat = 1;
		double[][] mspecdelta = null;
		double[][] mspecdeltadelta = null;

		System.out.println("Modulation Spectrogram matrix has size " + mspec.length + " X " + mspec[0].length);

		if (delta) {
			System.out.println("Adding deltas");
			mspecdelta = UtilsMath.computeDeltaMatrix(mspec);
			System.out.println("Added deltas " + mspecdelta.length + " X " + mspecdelta[0].length);
			nfeat++;
		}

		if (doubledelta) {
			System.out.println("Adding double deltas");
			mspecdeltadelta = UtilsMath.computeDoubleDeltaMatrix(mspec);
			System.out.println("Added double deltas " + mspecdeltadelta.length + " X " + mspecdeltadelta[0].length);
			nfeat++;
		}

		System.out.println("Building feature sequence");
		for (int i = 0; i < mspec.length; i++) {
			int j = 0;
			double[] featureRow = new double[nfeat * modulation_spectrogram_nfeatures];
			System.arraycopy(mspec[i], 0, featureRow, 0, modulation_spectrogram_nfeatures);
			j += modulation_spectrogram_nfeatures;
			if (delta) {
				System.arraycopy(mspecdelta[i], 0, featureRow, j, modulation_spectrogram_nfeatures);
				j += modulation_spectrogram_nfeatures;
			}
			if (doubledelta) {
				System.arraycopy(mspecdeltadelta[i], 0, featureRow, j, modulation_spectrogram_nfeatures);
				j += modulation_spectrogram_nfeatures;
			}
			features.add(featureRow);
		}

		System.out.println(
				"Feature list complete. Size " + mspec.length + " X " + (nfeat * modulation_spectrogram_nfeatures));

		return features;
	}

	public static double getTime(int index) {
		return ((double) index) * ModulationSpectrogram.windowShift;
	}

	public static void main(String[] args) throws Exception {
	    String usage =
	        "Usage: java ModulationSpectrogramManager -inputfile <input.wav> "
	        + "[-nfeats 8] [-maxfreq 3000] "
	        + "[-usedelta true] [-usedoubledelta true] "
	        + "[-outputfile <output.csv>]";

	    if (args.length == 0) {
	        System.err.println(usage);
	        return;
	    }

	    File signal = null;
	    File output = null;
	    int nFeatures = 8;
	    double maxFrequency = 3000.0;
	    boolean useDelta = true;
	    boolean useDoubleDelta = true;

	    for (int i = 0; i < args.length; i++) {
	        String option = args[i];

	        if ("-help".equals(option) || "--help".equals(option)) {
	            System.out.println(usage);
	            return;
	        }

	        switch (option) {
	            case "-inputfile":
	            case "-outputfile":
	            case "-nfeats":
	            case "-maxfreq":
	            case "-usedelta":
	            case "-usedoubledelta":
	                break;
	            default:
	                throw new IllegalArgumentException(
	                    "Unknown argument: " + option + "\n" + usage
	                );
	        }

	        if (i + 1 >= args.length || args[i + 1].startsWith("-")) {
	            throw new IllegalArgumentException(
	                "Missing value for " + option
	            );
	        }

	        String value = args[++i];

	        switch (option) {
	            case "-inputfile":
	                signal = new File(value);
	                break;

	            case "-outputfile":
	                output = new File(value);
	                break;

	            case "-nfeats":
	                nFeatures = Integer.parseInt(value);
	                break;

	            case "-maxfreq":
	                maxFrequency = Double.parseDouble(value);
	                break;

	            case "-usedelta":
	            case "-usedoubledelta":
	                if (!"true".equalsIgnoreCase(value)
	                        && !"false".equalsIgnoreCase(value)) {
	                    throw new IllegalArgumentException(
	                        option + " must be true or false."
	                    );
	                }

	                if ("-usedelta".equals(option)) {
	                    useDelta = Boolean.parseBoolean(value);
	                } else {
	                    useDoubleDelta = Boolean.parseBoolean(value);
	                }
	                break;
	        }
	    }

	    if (signal == null) {
	        throw new IllegalArgumentException(
	            "Required argument missing: -inputfile\n" + usage
	        );
	    }

	    if (!signal.isFile()) {
	        throw new IllegalArgumentException(
	            "Input file not found: " + signal.getAbsolutePath()
	        );
	    }

	    if (nFeatures <= 0) {
	        throw new IllegalArgumentException("-nfeats must be positive.");
	    }

	    if (!Double.isFinite(maxFrequency) || maxFrequency <= 0) {
	        throw new IllegalArgumentException(
	            "-maxfreq must be finite and positive."
	        );
	    }

	    if (output == null) {
	        output = new File(
	            signal.getAbsoluteFile().getParentFile(),
	            signal.getName().replaceFirst("(?i)\\.wav$", "")
	                + "_modulation_spectrogram.csv"
	        );
	    }

	    System.out.println(
	        "Calculating Modulation Spectrogram with " + nFeatures
	        + " features and a max frequency of " + maxFrequency
	    );
	    System.out.println("Use delta: " + useDelta);
	    System.out.println("Use double delta: " + useDoubleDelta);
	    System.out.println("Output file path: " + output.getAbsolutePath());

	    ModulationSpectrogramManager msm = new ModulationSpectrogramManager();

	    List<double[]> msFeatures = msm.extractModulationSpectrogram(
	        signal,
	        nFeatures,
	        maxFrequency,
	        useDelta,
	        useDoubleDelta
	    );

	    boolean header = true;
	    UtilsMath.saveFeaturesToFile(msFeatures, output, header);
	}
	
}
