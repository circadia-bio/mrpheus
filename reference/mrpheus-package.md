# mrpheus: Raw Physiological Signal Analysis for Biological Rhythms Research

Raw physiological signal analysis for biological rhythms research.
Provides ingestion, artefact detection, event detection, and feature
extraction for multi-modal physiological recordings including
polysomnography (EDF/EDF+), MRI-concurrent physiological logs (Philips
PMU), and EEG. Covers cardiac rhythm (QRS detection, HRV), respiratory
rhythm (apnoea detection, respiratory indices), and neural oscillations
(sleep spindles, slow oscillations, automatic AASM sleep staging via a
pre-trained LightGBM model ported from YASA; Vallat & Walker, 2021).
Exports staged hypnograms to hypnor and derived metrics to syncR. Part
of the Circadia Lab R ecosystem at Northumbria University.

## See also

Useful links:

- <https://mrpheus.circadia-lab.uk>

- <https://github.com/circadia-bio/mrpheus>

- Report bugs at <https://github.com/circadia-bio/mrpheus/issues>

## Author

**Maintainer**: Lucas França <lucas.franca@northumbria.ac.uk>

Authors:

- Lucas França <lucas.franca@northumbria.ac.uk>

- Mario Leocadio-Miguel ([ORCID](https://orcid.org/0000-0002-7248-3529))
