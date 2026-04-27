
<h1 align="center">RPA-wrapper</h1>
<a href="https://github.com/xph9876/RibosePreferenceAnalysis">Ribose Preferred Analysis</a> Wrapper Scripts

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>


[![Commits][Commits-shield]][Commits-url]
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Website][website-shield]][website-url]
[![Issues][issues-shield]][issues-url]
[![License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#recently-added">Recently Added</a></li>
    <li><a href="#installation">Installation</a></li>
      <ul>
        <li><a href="#getting-the-code">Getting the code</a></li>
        <li><a href="#creating-the-environment-with-required-dependencies">Creating the environment with required dependencies</a></li>
        <li><a href="#additional-dependencies">Additional Dependencies</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
      <ul>
        <li><a href="#defining-variables">Defining variables</a></li>
        <li><a href="#initializing-functions-and-activating-environment">Initializing functions and activating environment</a></li>
        <li><a href="#running-masterrunsh">Running Masterrun.sh</a></li>
        <li><a href="#output-structure">Output structure</a></li>
        <li><a href="#statistical-test-for-preference-and-genotype-comparisons">Statistical test for preference and genotype comparisons</a></li>
        <li><a href="#generating-stacked-barplots-for-hotspots">Generating Stacked barplots for hotspots</a></li>
      </ul>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#citations">Citations</a></li>
  </ol>
</details>

<!-- RECENTLY ADDED -->
## Recently Added

**Unified Masterrun.sh with parallel strand workflows**
`scripts/Masterrun.sh` now runs the both-strands, same-strand, and opposite-strand workflows simultaneously as background jobs, writing results into separate subdirectories (`both/`, `same/`, `opp/`) inside the configured output directory.

**Background frequency caching**
Background frequency computation (`bg_freq`, `bg_freq_ss`, `bg_freq_os`) is performed once and stored outside the per-run output directory. Subsequent runs with the same reference and range file skip this step automatically, saving significant compute time.

**Per-run output directories**
All output is written to a named directory set by the `outdir` variable in `Masterrun.sh`. Change `outdir` between runs to keep results from different datasets or parameters cleanly separated.

**Automatic timestamped logging**
Each workflow subshell logs all stdout and stderr to a `log/` directory inside its output subdirectory. Log files are named by timestamp (`YYYYMMDD_HHMMSS.log`) so reruns accumulate rather than overwrite. Background frequency logs are written to `<outdir>/log/`.

**Reorganized repository layout**
- `scripts/` — all shell and R scripts (`Masterrun.sh`, `Heatmapwrapper.sh`, `mww.R`, `mww_diff.R`, `comp.R`, `comp_bars.R`)
- `config/` — configuration files including the example `order` file
- `RPA/` — RibosePreferenceAnalysis Python scripts (git submodule, path resolved automatically)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- Installation -->
## Installation
### Getting the code
```sh
git clone --recurse-submodules https://github.com/twarner620/RPA-wrapper.git
```
### Creating the environment with required dependencies
```sh
conda env create --file environment.yml
```
### Additional Dependencies
* Input files (bed) containing single nucleotide locations, mainly for rNMP data. (another single nucleotide data can also be experimented on!)
* Reference genome files (.fa and .fai) of the organism being used (also used to generate bed files)
* Ranges BED file of the genome locations to be analyzed, for which background frequency will be calculated
* Order file — example provided in `config/order`

<!-- USAGE -->
## Usage
### Defining variables
Edit the top section of `scripts/Masterrun.sh`:
```bash
ref='path/to/reference/sacCer2/sacCer2.fa'  # Reference FASTA file
range='path/to/ranges/sacCer2/nucl.bed'      # Ranges BED file (nuclear genome)
#range='path/to/chrM.bed'
bed='path/to/bed/'                           # Folder of bed files
order='config/order'                         # Order file (example in config/)
outdir='run1'                                # Output directory; change per run to keep results separate
```
The RPA submodule path and all internal paths are resolved automatically — no `scripts` variable needed.

### Initializing functions and activating environment
```bash
conda activate RPAwrapper_env
source path/to/RPA-wrapper/scripts/Heatmapwrapper.sh
```

### Running Masterrun.sh
Edit the variables at the top of `scripts/Masterrun.sh`, then run from the repo root:
```bash
source path/to/RPA-wrapper/scripts/Masterrun.sh
```
This single command runs the full pipeline for all three strand modes in parallel:
- **bg_freq** (both, same, opposite strand) is computed once and cached in `bg_freq/` at the repo root. On subsequent runs with the same `ref` and `range`, this step is skipped.
- **Both-strand**, **same-strand**, and **opposite-strand** workflows then run simultaneously, each writing into its own subdirectory of `outdir`.

### Output structure
```
bg_freq/                        # Shared background frequency cache (ref+range specific)

<outdir>/
  log/                          # Background frequency logs (timestamped)
    bg_freq_YYYYMMDD_HHMMSS.log
    bg_freq_ss_YYYYMMDD_HHMMSS.log
    bg_freq_os_YYYYMMDD_HHMMSS.log
  both/                         # Both-strands results
    log/YYYYMMDD_HHMMSS.log
    sample_freq/
    norm_freq/
    plots/
  same/                         # Same-strand results
    log/YYYYMMDD_HHMMSS.log
    sample_freq/
    norm_freq/
    plots/
  opp/                          # Opposite-strand results
    log/YYYYMMDD_HHMMSS.log
    sample_freq/
    norm_freq/
    plots/
```

### Statistical test for preference and genotype comparisons
```bash
# Source Heatmapwrapper.sh and set variables as above, then from the relevant outdir subdir:
mww "$scripts" "$ref" "$range" "$bed" "$order"
```

### Generating Stacked barplots for hotspots
```bash
Rscript path/to/RPA-wrapper/scripts/comp.R -m sorted_chrM_mono_0 #hotspot composition files usually contain one entry for every genotype.
```


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GNU GPL3 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact
Deepali L. Kundnani - [deepali.kundnani@gmail.com](mailto::deepali.kundnani@gmail.com)    [![LinkedIn][linkedin-shield]][linkedin-url] 
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Citations
Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!+
* <b> Distinct features of ribonucleotides within genomic DNA in Aicardi-Goutières syndrome (AGS)-ortholog mutants of <i>Saccharomyces cerevisiae</i> </b>
Deepali L. Kundnani, Taehwan Yang, Alli L. Gombolay, Kuntal Mukherjee, Gary Newnam, Chance Meers, Zeel H. Mehta, Celine Mouawad, Francesca Storici
bioRxiv 2023.10.02.560505; doi:[https://doi.org/10.1101/2023.10.02.560505]( https://doi.org/10.1101/2023.10.02.560505)
* Kundnani, D. (2024). rNMP_hotspots:2.0.0 (2.0.0). Zenodo.  [https://doi.org/10.5281/zenodo.8152090](https://doi.org/10.5281/zenodo.8152090) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8152090.svg)](https://doi.org/10.5281/zenodo.8152090)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/DKundnani/RPA-wrapper?style=for-the-badge
[contributors-url]: https://github.com/DKundnani/RPA-wrapper/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/DKundnani/RPA-wrapper?style=for-the-badge
[forks-url]: https://github.com/DKundnani/RPA-wrapper/forks
[stars-shield]: https://img.shields.io/github/stars/DKundnani/RPA-wrapper?style=for-the-badge
[stars-url]: https://github.com/DKundnani/RPA-wrapper/stargazers
[issues-shield]: https://img.shields.io/github/issues/DKundnani/RPA-wrapper?style=for-the-badge
[issues-url]: https://github.com/DKundnani/RPA-wrapper/issues
[license-shield]: https://img.shields.io/github/license/DKundnani/RPA-wrapper?style=for-the-badge
[license-url]: https://github.com/DKundnani/RPA-wrapper/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/deepalik
[product-screenshot]: images/screenshot.png
[commits-url]: https://github.com/DKundnani/RPA-wrapper/pulse
[commits-shield]: https://img.shields.io/github/commit-activity/t/DKundnani/RPA-wrapper?style=for-the-badge
[website-shield]: https://img.shields.io/website?url=http%3A%2F%2Fdkundnani.bio%2F&style=for-the-badge
[website-url]:http://dkundnani.bio/ 
