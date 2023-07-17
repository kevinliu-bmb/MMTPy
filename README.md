# MMTpy - Microbiome Modeling Toolbox Extensions

<!-- TOP OF README.MD -->
<a name="readme-top"></a>

<!-- PROJECT SHIELDS -->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Issues][issues-shield]][issues-url]
[![GNU v3.0 License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <h3 align="center">MMTpy</h3>

  <p align="center">
    Microbiome Modeling Toolbox extensions and workflows for multi-species genome-scale models of metabolism
    <br />
    <a href="https://github.com/kevinliu-bmb/MMTpy/issues">Report Bug</a>
    ·
    <a href="https://github.com/kevinliu-bmb/MMTpy/issues">Request Feature</a>
    ·
    <a href="https://zomorrodi.mgh.harvard.edu">The Zomorrodi Lab @ MGH</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details open>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#a-note-regarding-model-compartments-and-reaction-types">A note regarding model compartments and reaction types</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

Python-based tools for the processing and analysis of multi-species GEnome-scale Models (GEM) of Metabolism that is compatible with multi-species community models constructed using the [mgPipe.m pipeline of the Microbiome Modeling Toolbox](https://opencobra.github.io/cobratoolbox/latest/modules/analysis/multiSpecies/microbiomeModelingToolbox/index.html) by [Heinken et al. (2022)](https://academic.oup.com/bioinformatics/article/38/8/2367/6528309).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

The following information will guide you through installing the required Python module prerequisites.
To get a local copy up and running, follow these simple example steps.

### Prerequisites

_Several core functionalities of MMTpy rely on the [COBRApy](https://github.com/opencobra/cobrapy) and [OptLang](https://github.com/opencobra/optlang) Python modules; please refer to ```requirements.txt``` for the full list of dependencies and/or use the following guide to install all required dependencies._

1. Clone the repo

   ```sh
   git clone https://github.com/kevinliu-bmb/MMTpy.git
   ```

2. Install required Python modules (the project was created under a conda environment using Python 3.9.16)

   ```sh
   pip install -r requirements.txt
   ```

***Important Note: for improved computational efficiency, please install [Gurobi](https://www.gurobi.com) or other compatible LP solvers during optimization tasks, such as the [IBM CPLEX Optimizer](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer), and obtain a valid license. Free academic licenses for Gurobi can be obtained through the [Academic Program and Licenses](https://www.gurobi.com/academia/academic-program-and-licenses/) page.***

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

To run a single instance of the FBA simulation workflows using COBRApy, import the ```mmtpy_workflows.py``` script using a Python console (tested on Python 3.9.16) under the cloned GitHub repository folder, as shown in the two examples below.

   ```python
   > from mmtpy_opt_workflows import optimize_model, optimize_model_mbx
   > # run naive FBA simulations on a model already loaded into memory.
   > optimize_model(model_input="example_data/models/microbiota_model_diet_Case_1_18_month.json", output_path="example_outputs")
   > # model input can also be a path to a .mat or .json file.
   > optimize_model(model_input=model, output_path="example_outputs")
   > # run metabolomics data-constrained FBA simulations.
   > optimize_model_mbx(model_input=model, mbx_path="example_data/metabolomics_data.csv", output_path="example_outputs")
   ```

To launch multiple instances of both FBA simulation workflows simultaneously for more than one model in parallel, configure the relevant paths in the ```run_workflows_parallel.py``` script, save the script, and run it in the command line, as shown in the example below. *This option is not recommended for models with large numbers of reactions and/or metabolites, as the 'optimize_model' workflow tends to be more computationally intensive than the 'optimize_model_mbx' workflow and tends to result in a longer runtime.*

   ```sh
   python run_workflows_parallel.py
   ```

To launch multiple instances of a single FBA simulation workflow simultaneously for more than one model in parallel, configure the relevant paths and workflow type in the ```run_single_workflows_parallel.py``` script, save the script, and run it in the command line, as shown in the example below. *This option is not recommended for runs that contain a large number of input models, as the script will attempt to load all models into memory simultaneously, which may result in memory errors.*

   ```sh
   python run_single_workflows_parallel.py
   ```

Other convenient tools, such as ```set_default_bounds()``` for resetting the model reactions bounds and ```convert_model_format()``` to convert any COBRApy-supported model format to JSON format, can be called within a Python console after importing ```cobra_utils.py```.

   ```python
   > from mmtpy_utils import set_default_bounds, convert_model_format
   > # set the model bounds to MMTpy default conventions.
   > set_default_bounds(model=model, source="MMTpy")
   > # convert a COBRApy supported model format to JSON format.
   > convert_model_format(model_path=model, output_path="example_data/models")
   ```

An additional feature available in MMTpy is the ability to match metabolite names from GC-MS (or alternative instrument) quantified and annotated output metabolite names, which can include common names or any other non-standard biochemical nomenclature, to VMH metabolite identifiers through the ```cobra_utils.py``` script and the included exhaustive list of VMH metabolites and their respective alternative identifiers in ```all_vmh_metabolites.tsv```. The usage of the VMH metabolite identifier database enables several metabolite matching strategies, such as through InChIString, InChIKey, CID, and isomeric SMILES (only used as a last resort due to the possibility of stereoisomers), found under the ```~/data_dependencies/``` directory. As a fallback strategy, a manually curated mapping file is also provided as ```manually_matched_keys.txt```, which enables the usage of the mapping function in the absence of internet access in addition to providing a more comprehensive mapping of GC-MS names to VMH identifiers. The ```match_names_to_vmh``` function can be called within a Python console after importing ```cobra_utils.py```, as shown in the example below.

   ```python
   > from mmtpy_utils import *
   > match_names_to_vmh(model_input=model, mbx_filepath="example_data/metabolomics_data.csv", output_filepath="example_outputs")
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- A NOTE REGARDING MODEL COMPARTMENTS AND REACTIONS -->
## A Note Regarding Model Compartments and Reaction Types

***The following information is based on inferred definitions from [Microbiome Modeling Toolbox](https://opencobra.github.io/cobratoolbox/latest/modules/analysis/multiSpecies/microbiomeModelingToolbox/index.html) and may not represent the original author's definitions. Furthermore, these definitions are unlikely to be generalizable to alternative tools. Nonetheless, we provide the following details if one were to apply MMTpy functions to models constructed using other tools.***

### Types of Compartments

#### There are five types of compartments used in GEMs

- [c] - Cytoplasmic space
- [p] - Periplasmic space
- [u] - Extracellular space (often denoted as [e] in other models)
- [fe] - Fecal
- [d] - Diet

##### Compartment Details

- The [c] compartment is the cytoplasmic space, where most of the microbe's internal metabolic reactions occur.
- The [p] compartment is the periplasmic space, which is the space between the inner and outer membranes of Gram-negative bacteria.
- The [u] compartment is the extracellular space, which is the space outside of the cell.
- The [fe] compartment is the fecal compartment, which is where the fecal metabolites are present.
- The [d] compartment is the diet compartment, which is where the diet metabolites are present.

### Types of Reactions

#### Diet exchange reactions

- 'Diet_EX_met[d]': 'met[d] <=>’
- Reversible reaction
- Typically only has negative bounds with values that are defined by the diet.

#### Diet transport reactions (i.e., "DUt" reactions)

- 'DUt_met': 'met[d] --> met[u]'
- Irreversible reaction
- By default, has bounds of (0., 1000000.)

#### Internal exchange reactions (i.e., "IEX" reactions)

- 'microbe_IEX_met[u]tr': 'microbe_met[c] <=> met[u]'
- Reversible reaction
- By default, has bounds of (-1000., 1000.)

#### Fecal transport reactions (i.e., "UFEt" reactions)

- 'UFEt_met': 'met[u] --> met[fe]'
- Irreversible reaction
- By default, has bounds of (0., 1000000.)

#### Fecal exchange reactions

- 'EX_met[fe]': 'met[fe] <=>'
- Reversible reaction
- Has bounds of (-1000., 1000000.) according to the [Microbiome Modeling Toolbox](https://opencobra.github.io/cobratoolbox/latest/modules/analysis/multiSpecies/microbiomeModelingToolbox/index.html)* and bounds of (0., 1000000.) according to MMTpy conventions

  *Note: fecal exchange reactions for the "microbeBiomass" have bounds of (-10000., 1000000.) according to the [Microbiome Modeling Toolbox](https://opencobra.github.io/cobratoolbox/latest/modules/analysis/multiSpecies/microbiomeModelingToolbox/index.html) and is defaulted to (0., 1000000.) in MMTpy.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

- [x] Add FBA simulations feature using the [COBRApy](https://github.com/opencobra/cobrapy) module
- [x] Add non-standardized metabolite name mapping to [Virtual Metabolic Human (VMH)](https://www.vmh.life/) identifiers feature
- [x] Add metabolomic data-constrained optimization feature using the [COBRApy](https://github.com/opencobra/cobrapy) module
- [ ] Add multi-species model creation with gap-filling feature

See the [open issues](https://github.com/kevinliu-bmb/MMTpy/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/yourFeature`)
3. Commit your Changes (`git commit -m 'Add some yourFeature'`)
4. Push to the Branch (`git push origin feature/yourFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the GNU v3.0 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Kevin Liu - [GitHub](https://github.com/kevinliu-bmb) - <kliu@hms.harvard.edu> - Department of Biomedical Informatics, Harvard Medical School

Project Link: [https://github.com/kevinliu-bmb/MMTpy](https://github.com/kevinliu-bmb/MMTpy)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

MMTpy was created based on the Microbiome Modeling Toolbox by [Almut Heinken](https://scholar.google.com/citations?user=4Lu-c34AAAAJ&hl=en&oi=ao) and [Ines Thiele](https://orcid.org/0000-0002-8071-7110).

The tools here were conceptualized by [Ali R. Zomorrodi, PhD.](https://orcid.org/0000-0002-9134-8082), and its elementary implementation was done by [Izzy Goodchild-Michelman](https://www.linkedin.com/in/isabella-goodchild-michelman-921bab196).

MMTpy was developed by [Kevin Liu](https://dbmi.hms.harvard.edu/people/kevin-liu) within the [Zomorrodi lab](https://zomorrodi.mgh.harvard.edu/) at Massachusetts General Hospital under the mentorship of the Principal Investigator, [Ali R. Zomorrodi, PhD](https://orcid.org/0000-0002-9134-8082).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
[contributors-shield]: https://img.shields.io/github/contributors/kevinliu-bmb/MMTpy.svg?style=for-the-badge
[contributors-url]: https://github.com/kevinliu-bmb/MMTpy/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/kevinliu-bmb/MMTpy.svg?style=for-the-badge
[forks-url]: https://github.com/kevinliu-bmb/MMTpy/network/members
[issues-shield]: https://img.shields.io/github/issues/kevinliu-bmb/MMTpy.svg?style=for-the-badge
[issues-url]: https://github.com/kevinliu-bmb/MMTpy/issues
[license-shield]: https://img.shields.io/github/license/kevinliu-bmb/MMTpy.svg?style=for-the-badge
[license-url]: https://github.com/kevinliu-bmb/MMTpy/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/kevin-liu-
