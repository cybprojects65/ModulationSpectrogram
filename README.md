# Modulation Spectrogram
An implementation of the Modulation Spectrogram feature extraction process from Greenberg and Kingsury (1997).

Main class: **it.cnr.speech.filters.ModulationSpectrogram**


	it.cnr.speech.modulationspec.main.ModulationSpectrogramManager -nfeats 8 -maxfreq 3000 -usedelta true -usedoubledelta true -inputfile "samples/PS14Audio_noisy_speech.wav" -outputfile "samples/PS14Audio_noisy_speech_modspec.csv"
	
	java -cp ./ -jar modulation_spectrogram.jar -nfeats 8 -maxfreq 3000 -usedelta true -usedoubledelta true -inputfile "./samples/PS14Audio_noisy_speech.wav" -outputfile "./samples/PS14Audio_noisy_speech_modspec.csv"
