# Modulation Spectrogram

A Java implementation of the modulation spectrogram feature extraction
process described by Greenberg and Kingsbury (1997), “The modulation
spectrogram: In pursuit of an invariant representation of speech.”

The implementation extracts features describing slow amplitude variations
within acoustic frequency bands, focusing on modulation around 4 Hz.

## Entry points

- **Command-line entry point:**  
  `it.cnr.speech.modulationspec.main.ModulationSpectrogramManager`
- **Core feature extraction:**  
  `it.cnr.speech.filters.ModulationSpectrogram`

## Usage

### Eclipse

Set the main class to:

it.cnr.speech.modulationspec.main.ModulationSpectrogramManager

Set the program arguments to:

-nfeats 8 -maxfreq 3000 -usedelta true -usedoubledelta true -inputfile "samples/PS14Audio_noisy_speech.wav" -outputfile "samples/PS14Audio_noisy_speech_modspec.csv"

### Executable JAR

java -jar modulation_spectrogram.jar -nfeats 8 -maxfreq 3000 -usedelta true -usedoubledelta true -inputfile "./samples/PS14Audio_noisy_speech.wav" -outputfile "./samples/PS14Audio_noisy_speech_modspec.csv"

### Arguments

| Argument | Description | Default |
|---|---|---|
| `-inputfile` | Input audio file | Required |
| `-outputfile` | Output CSV file | Input basename followed by `_modulation_spectrogram.csv` |
| `-nfeats` | Number of acoustic frequency bands and base features | `8` |
| `-maxfreq` | Upper acoustic filter-bank boundary in Hz | `3000` |
| `-usedelta` | Include first-order temporal delta features | `true` |
| `-usedoubledelta` | Include second-order temporal delta features | `true` |
| `-help` or `--help` | Display usage information | — |

`-maxfreq` controls the acoustic frequency range, not the modulation
frequency. Modulation features are evaluated at 4 Hz.

## Supported input

Use signed, 16-bit, little-endian, mono PCM WAV files.

A sampling rate of **16 kHz** is suitable. The processing pipeline also
accepts other whole-number sampling rates of at least **8 kHz**, subject
to the supported PCM format, and resamples the audio to 8 kHz internally.
Input already sampled at 8 kHz is left at that rate.

Stereo and other sample encodings must be converted before processing.

## Processing pipeline

For each recording, the implementation:

1. Loads the audio and resamples it to 8 kHz, applying anti-aliasing
   filtering when reducing the sampling rate.
2. Separates the signal into acoustic frequency bands.
3. Applies half-wave rectification to each band.
4. Applies a 28 Hz low-pass filter to extract a slowly varying envelope.
5. Downsamples each envelope to 80 Hz.
6. Normalizes each envelope by its average absolute level over the recording.
7. Applies 20-sample Hamming windows, corresponding to 250 ms, with a
   one-sample hop, corresponding to 12.5 ms.
8. Directly calculates the 4 Hz DFT coefficient and squares its magnitude.
9. Converts the power to dB and clips it to the configured implementation
   range of −30 to +30 dB in the manager's extraction path.
10. Optionally calculates delta and double-delta features and exports CSV.

Only complete modulation-analysis windows are included.

## Code structure

### `it.cnr.speech.modulationspec.main`

**`ModulationSpectrogramManager`**

Provides the command-line interface and coordinates feature extraction.
It parses arguments, invokes the core processor, arranges results as
frames × features, appends optional derivatives, and writes the CSV.

Its `extractModulationSpectrogram(...)` method can also be called directly
from another Java application.

### `it.cnr.speech.filters`

**`ModulationSpectrogram`**

Coordinates the processing stages for each acoustic band. Internally,
the base modulation spectrogram is stored as bands × frames.

The `saveSteps` flag enables intermediate WAV exports for inspection.

**`GreenwoodFilterBank`**

Defines the acoustic band weights and applies them through the inherited
FFT-based filtering machinery. Despite the class name, the current
frequency-spacing calculation uses the Mel scale.

**`LowPassFilterDynamic`**

Provides windowing, FFT, inverse FFT, and overlap-add operations used by
the acoustic filter bank. The revised 28 Hz envelope filtering is performed
separately by `SignalProcessing.lowPassEnvelope(...)`.

**`HalfWaveRectification`**

Replaces negative signal samples with zero before envelope filtering.

**`DownSampler`**

Performs integer-factor decimation by retaining regularly spaced samples.
Anti-aliasing filtering must be performed before calling it.

**`Delta` and `HighPassFilter`**

Additional utility classes, not used by the manager's current extraction
pipeline.

### `it.cnr.speech.utils`

**`SignalProcessing`**

Provides audio loading, initial resampling, envelope normalization,
Butterworth envelope filtering, and direct calculation of the 4 Hz
modulation coefficient. It also retains general FFT-related helpers.

**`UtilsMath`**

Provides matrix transposition, temporal delta and double-delta
calculations, and CSV export with frame timestamps.

**`AudioBits`**

Reads audio samples and exposes their format information.

**`AudioWaveGenerator`**

Writes intermediate audio signals to WAV files.

## Output format

The manager writes one CSV row per analysis frame:

Time_start_s,Time_end_s,F0,F1,...

The feature columns contain, in order:

1. Base modulation features.
2. Delta features, when enabled.
3. Double-delta features, when enabled.

Base bands are currently exported from higher to lower acoustic frequency.
The derivative blocks preserve that ordering.

With eight bands and both derivative options enabled, the output contains
**24 feature columns plus two timestamp columns**.

Frame starts are spaced 12.5 ms apart, and each frame covers a 250 ms
analysis window. The causal envelope filter introduces delay; timestamps
describe the analysis windows and do not compensate for this delay.

## Implementation notes

This implementation follows the paper's processing concept but is not
an exact reproduction. In particular, it uses Mel-spaced acoustic bands,
a 16th-order Butterworth envelope filter, and fixed dB clipping rather
than the paper's peak-relative display scaling.

The 20-sample window, one-sample hop, and 4 Hz coefficient are currently
fixed in `modulationMagnitude4Hz(...)`. Changing the older window-related
fields alone does not change that calculation.

## Build dependencies

- Java 16 or later.
- Apache Commons Math 3.6.1.

The Maven configuration uses `src` as the source directory. The repository
also includes an executable JAR and an Eclipse-generated Ant export script.
The Ant script contains machine-specific paths that must be adapted before
use on another system.
	