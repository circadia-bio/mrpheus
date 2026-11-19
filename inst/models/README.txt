Pre-trained YASA sleep staging classifiers (Vallat & Walker, 2021), converted to
gzip-compressed LightGBM text files by data-raw/fetch_yasa_models.py:

  yasa_eeg.txt.gz          yasa_eeg_demo.txt.gz
  yasa_eeg_eog.txt.gz      yasa_eeg_eog_demo.txt.gz
  yasa_eeg_emg.txt.gz      yasa_eeg_emg_demo.txt.gz
  yasa_eeg_eog_emg.txt.gz  yasa_eeg_eog_emg_demo.txt.gz

stage_epochs() chooses the file that matches the channels it is given (and
age / male for the "demo" files), as yasa.SleepStaging does. Source: YASA
classifiers clf_*_lgb_0.5.0.joblib (BSD 3-Clause; see README.md).
