# IsoQ - Single-sample metabolomics and isotopic analysis

## What is IsoQ?
**IsoQ is a scientific software dedicated to the automated processing of metabolomics and isotopic MS measurements**.
IsoQ corrects raw MS data (mass fractions) for
naturally-occurring isotopes of all elements and purity of the
isotopic tracer using IsoCor, and calculates metabolite concentrations 
from calibration curves.

It is one of the routine tools that we use at the [MetaSys team](http://www.lisbp.fr/en/research/molecular-physiology-and-metabolism/metasys.html) and [MetaToul platform](http://www.metatoul.fr) in isotopic studies of metabolic systems.

Additional information can be found in the following [publication](https://doi.org/SUBMITTED):

*Simultaneous measurement of absolute metabolite concentration and isotope incorporation by mass spectrometry*, by Maud Heuillet, Pierre Millard, Madi Cissé, Matthieu Guionnet, Laetitia Linares, Fabien Létisse, Laurent Le Cam, Jean-Charles Portais and Floriant Bellvert.


The code is open-source, and available under a GPLv3 license.

## Installation
IsoQ requires Python 3.5 or higher and run on all plate-forms. If you do not have a Python environment
configured on your computer, we recommend that you follow the instructions
from [Anaconda](https://www.anaconda.com/download/).

To install IsoQ locally, open a terminal (e.g. run *Anaconda Prompt* if you have installed Anaconda) and type:

```bash
$ pip install git+https://github.com/MetaSys-LISBP/IsoQ.git
```

Alternatively, you can download the github repository and run:

```bash
$ pip install /path/to/IsoQ/
```

An example of IsoQ usage is provided below and in `*/isoq/data/example_process.py*`.

## Usage

First, load IsoQ in your Python session (IsoQ must be installed beforehand, as detailed above):

```python
> import isoq.process as iq
```

The processing of a dataset starts by defining some information relative to the dataset (e.g. data directory) as well as the processing parameters:

Data files
```python
> datafile = "example_dataset.csv"
> calibfile = "example_calibration.csv"
> data_folder = "C:/Users/millard/Documents/GIT/IsoQ/IsoQ/isoq/data/"
```

Processing parameters
```python
> # isotopic tracer
> tracer = '13C'
> # resolution of the MS instrument, m/z at which resolution is measured, and type of instrument (see IsoCor documentation for details)
> resolution = 70000
> mz_of_resolution = 400
> resolution_formula_code = 'orbitrap'
> # options relative to the purity of the isotopic tracer (see IsoCor documentation for details)
> tracer_purity = [0.01, 0.99]
> correct_NA_tracer = True
> # isotopic purity of the other isotope used in the IS
> purity15N = [0.01, 0.99]
```

Then, process your data by running the following command:

```python
> iq.run(data_folder, # data directory (str)
>        datafile, # name of the file containing the samples data (str)
>        calibfile, # name of the file containing the calibration data (str)
>        tracer, # isotopic tracer (str)
>        resolution, # resolution of the MS spectrometer (float)
>        mz_of_resolution, # m/z at which the resolution is given (float)
>        tracer_purity, # isotopic purity of the tracer (vector)
>        correct_NA_tracer, # correct (or not) for natural abundance of the tracer element (boolean)
>        resolution_formula_code, # formula to calculate resolution on the MS instrument (str)
>        purity15N=purity15N, # 15N purity of the 13C15N-IS (vector)
>        verbose=False # verbose output (boolean)
>        )
```

The following information are saved in the data directory:

- `results.pdf`: contains the calibration results (calibration curves, R², relative residuals, etc)
- `results_CID.csv`: contains the isotopic data (carbon isotopologue distributions, mean enrichment, residuum, etc) as detailed in IsoCor documentation
- `results_conc.csv`: contains the calculated metabolite concentrations
- `results.log`: contains detailed information about the processing (correction parameters, warnings and errors that may occur during processing, etc)

## Bug and feature requests
If you have an idea on how we could improve IsoQ please submit a new *issue*
to [our GitHub issue tracker](https://github.com/MetaSys-LISBP/IsoQ/issues).

## Contributions
Contributions are very welcome! :heart:

In development mode, do a `pip install -e /path/to/IsoQ/` to install locally the development version.

Please work on your own fork and
follow [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide.

## How to cite


Heuillet M., Millard P., Cissé M., Guionnet M., Linares L., Létisse F., Le Cam L., Portais J.-C. and Bellvert F. Simultaneous measurement of absolute metabolite concentration and isotope incorporation by mass spectrometry. Submitted.


## Authors
Pierre Millard, Baudoin Delépine 

## Contact
:email: Pierre Millard, millard@insa-toulouse.fr
