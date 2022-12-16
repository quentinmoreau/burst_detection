# Burst Detection Algorithm - usage notes

library accompanying the preprint "Diverse beta burst waveform motifs characterize movement-related cortical dynamics" 

<https://www.biorxiv.org/content/10.1101/2022.12.13.520225v1>

## PYTHON

### DEPENDENCIES

MNE-Python (only for filtering)

scipy

numpy

### INSTALLATION

1. install dependencies with `pip` or `conda install`

2. add folder to the PYTHONPATH or put it in the same directory as the script


### USAGE

```python
from burst_detection import extract_bursts
```




```python
bursts = extract_bursts(
    raw_trials, tf, times, search_freqs, 
    band_lims, aperiodic_spectrum, sfreq, w_size=.26
)
```




```python
bursts_single_trial = extract_bursts_single_trial(
    raw_trial, tf, times, search_freqs, 
    band_lims, aperiodic_spectrum, sfreq, w_size=.26
)
```

>A single trial burst extraction. Single trial time-course, and single trial 
time-frequency. No regressing out the average ERF. Potentially 7-10 Hz high-pass 
of the time domain trial can get rid of evoked response related burst shape distortions.



#### Arguments:

`raw_trials: raw data for each trial (trial x time)`

>Array of trials from a single channel in time domain.


\
`tf: time-frequency decomposition for each trial (trial x freq x time)`

>Array of time-frequency data from a single channel. If the targeted range 
of frequencies is 13-30 Hz, it is recommended to select the TF with at least 3 
Hz buffer (10-33 Hz). Recommended use of superlet (Moca et al., 2021) 
transformed data.

\
`times: time steps`

>Time points of the `raw_trials` and `tf` should be equivalent


\
`search_freqs: frequency limits to search within for bursts (should be wider than band_lims)`

>List of frequencies corresponding to the selected ones in `tf`.


\
`band_lims: keep bursts whose peak frequency falls within these limits`

>The actual target frequency range. Bursts with peak frequency beyond this range
will be discarded.


\
`aperiodic_spectrum: aperiodic spectrum`

>Assuming the PSD was calculated based on the TF data averaged over time, the
frequency resolution of the PSD and TF should be the same. 


\
`sfreq: sampling rate`

>Assuming the sampling rate is the same for `raw_trials` and `TF`


\
`w_size=.2: size of the burst window in time domain`

>Default argument. Window size based on lagged coherence of the MEG data.



#### Output:

The function returns a dictionary with time-frequency burst features, and time 
domain waveforms. 

```python
{
    'trial': [],
    'waveform': [],
    'peak_freq': [],
    'peak_amp_iter': [],
    'peak_amp_base': [],
    'peak_time': [],
    'peak_adjustment': [],
    'fwhm_freq': [],
    'fwhm_time': [],
    'polarity': [],
    'waveform_times': []
}
```


`trial` index of a trial where the burst was detected

`waveform` array of burst waveforms [burst x time]

`peak_freq` peak frequency of the burst

`peak_amp_iter` relative TF amplitude of the burst during the iterations

`peak_amp_base` absolute TF amplitude of the burst (with the aperiodic
 spectrum subtracted)

`peak_time` TF peak time of the burst

`peak_adjustment` adjustment of the peak in ms

`fwhm_freq` frequency span of the burst

`fwhm_time` duration of the burst

`polarity` 0 - the polarity was not flipped, 1 - polarity was flipped

`waveform_times` 1d array containing the timepoints for a waveform


## MATLAB

### DEPENDENCIES

Fieldtrip


### USAGE
```
bursts = extract_bursts(raw_trials, tf, times, search_freqs, band_lims, aperiodic_spectrum, sfreq)
```

```
bursts = extract_bursts_single_trial(raw_trial, tf, times, search_freqs, band_lims, aperiodic_spectrum, sfreq)
```

#### Arguments:

AS ABOVE

#### Output:

AS ABOVE
