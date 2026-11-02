# PFSimple

## General information

The Paticle Finder Simple package is simplified version of the KF Partice package based on its mathematical apparatus. It is developed for the complete reconstruction of short-lived particles with their momentum, energy, mass, lifetime, decaylength, rapidity, etc. 2- and 3-body decays can get reconstructed. It is also possible to reconstruct multiple decays in one run as well as cascade decays.

## Requirements

### Root

https://root.cern/install/

Version 6.20 or newer is recommended. Root with c++17 standard is strongly recommended.

### AnalysisTree

Optional. Needed for building AnalysisTree interface.

https://github.com/HeavyIonAnalysis/AnalysisTree

Both pre-installed AnalysisTree and automatically installed together with PFSimple framework can be used (see chapter Installation)

AnalysisTree with c++17 standard is strongly recommended.

## Installation

Clone PFSimple
    
Source ROOT
```
source /path-to-root/install/bin/thisroot.sh
```
To use pre-installed AnalysisTree set the environment variable:
```
export AnalysisTree_DIR=/path-to-analysistree/install/lib/cmake/AnalysisTree
```
To build AnalysisTree automatically together with PFSimple use following cmake keys:
```
-DPFSimple_BUNDLED_AT=ON
-DPFSimple_BUNDLED_AT_VERSION=v2.3.4
```

To patch the bundled KFParticle so that it is built **without** `-ffast-math`
(independent of the chosen build type), add:
```
-DPFSimple_KFParticle_NO_FAST_MATH=ON
```
This is **OFF by default**, so by default KFParticle is built exactly as upstream
(with whatever flags its own CMakeLists.txt picks for the build type, including
`-ffast-math` for a `Release` build). It is only needed when PFSimple must
reproduce the numerical results of a reference build that was compiled without
`-ffast-math` (e.g. when cross-checking against the `afkf` package). The patch is
applied to the KFParticle source at configure time, so enable it once and build
from scratch (or remove the KFParticle sub-build) for it to take effect.
    
You need to source root and export AnalysisTree each time when you are compiling project from 0 (perform cmake command) but have no need to do it when just recompiling project (perform just make).
    
Install PFSimple
```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-pfsimple /path-to-source-pfsimple
make -j install
```
## Configuration of decay settings
The program is configured with a JSON config file which contains various settings related to input and output as well as the decays to reconstruct. In the following, the different blocks and fields are explained. In `at_interface/config/config_template.json` all available JSON fields are listed.

### io
| Key                         | Description                                                                                                                                                                   |
| ----------------------------| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input_treename`            | name of the input AnalysisTree (for standard analysistree `rTree`, after running PID framework `pTree`)                                                                       |
| `rectracks_branchname`      | branchname of reconstructed tracks in input analysistree, default names: `VtxTracks` for standard analyistree, `RecParticles` after running Pid-framework                     |
| `receventheader_branchname` | branchname of reconstructed event header in input analysistree, default names: `RecEventHeader` for standard analyistree, `AnaEventHeader` after running centrality-framework |
| `n_events`                  | number of events to be processed, set to `-1` to process all events                                                                                                           |
| `mc_info_available`         | Indicates if the data comes with mc info and in case it does, computes additional MC-based fields (set to `false` when analyzing Verticla Test data)                          |
| `mc_match_weight_min` (opt.)| Optional minimum MC-match weight for `pid_mode 1` (default: `0` = no threshold). When set `> 0`, a track is dropped unless it has an MC match whose `match_weight` field is `>= mc_match_weight_min`. Mirrors the `afkf` track selection, which drops tracks whose MC match weight is below `0.7`. Off by default, so standard `pid_mode 1` behaviour is unchanged. |
| `save_options` (optional)   | defines save options for the candidates of all decays in this config (see next table)                                                                                         |

| `save_options` Flag         | Description                                                                                                                                                                   |
| ----------------------------| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `make_plain_tree`           | normal root tree containing the candidates is written besides the default output tree in AnalysisTree format                                                                  |
| `signal_only`               | save only mc-true signal candidates                                                                                                                                           |
| `write_detailed_bg`         | detailed background information in output variable `generation`*                                                                                                              |

*Format of detailed background information in the `generation` variable:

format for bg for 2-body-decay: -m12 0 d2 d1  
format for bg for 3-body-decay: -m13 m23 m12 d3 d2 d1  
with m12: mother of daughter 1 & 2 etc.  
code for daughters (d):  
1 - reco daughter is unmatched to mc  
2 - reco daughter is matched, but primary  
3 - reco daughter is secondary, produced not in decay from mother with not expected pdg  
4 - reco daughter is secondary, produced not in decay from mother with expected pdg  
5 - reco daughter is secondary, produced in decay from mother with not expected pdg  
6 - reco daughter is secondary, produced in decay from mother with expected pdg  
code for mothers (m):  
1 - daughters have the same mother  
2 - daughters have different mothers  
0 - at least one daughter does not have mother (e.g. primary)  

If detailed background is not selected or if mc match for mother is found, value of variable "generation" is:  
0 - no mc match for mother  
1 - mother is primary particle  
2 - mother is second generation etc.

### pid_mode
Methods for particle identification:
| `pid_mode` | Description                                                                                               |
| ---------- | --------------------------------------------------------------------------------------------------------- |
| `0`        | no pid (only charge information is used! optional daughers must be set to 1 / -1 for pos / neg pdgs!)     |
| `1`        | mc-pid                                                                                                    |
| `2`        | reconstructed TOF pid - default from Pid-framework*                                                       |
| `3`        | reconstructed TOF pid - pdg with max. pdg-purity is selected, if pdg-purity >min. purity*                 |
| `4`        | reconstructed TOF pid - pdg is selected, if pdg-purity > min. purity (pdg-specfic purities are possible)* |

** *for 2-4: Pid-framework needs to be applied first to include TOF-pid & probabilities into the analysistree input-file**

### pid_purity
Mimimum purities for pid-mode 3 & 4:
- `default`:  minimum purity is set for all pdgs (pdg-spefic purities need to be omitted) (in pid-mode 3 & 4)
- `protons`,  `pions`, etc:  minimum purity is set specifically for every pdg (general purity will be overwritten, if pdg-specific purity is given) (in pid-mode 4)

### accept_fragments
Optional top-level flag (default: `false`). When set to `true`, nuclear fragments (`|pdg| >= 1e9`, e.g. He3/He4) are accepted as proton-hypothesis positive daughters, mirroring `afkf`'s `Find2DaughterDecay`. Fragments are always secondary, so they bypass the `chi2prim >= cut` gate. This is only needed to reproduce `afkf`'s fragment acceptance in a one-off `afkf`-vs-PFSimple comparison; leave it off (`false`, the default) otherwise.

### chi2_prim_afkf
Optional top-level flag (default: `false`). When set to `true`, the chi2-to-primary-vertex computation is done the way the `afkf` package does it (scalar `KFParticle(ToKFPTrack(track), pdg).GetDeviationFromVertex(pv)`), instead of the default legacy SIMD path. The legacy SIMD path adds an extra `TransportToPoint` that lowers `chi2_prim` and undercounts secondary tracks. This is only needed to reproduce `afkf`'s `ComputeChiToVertex` exactly in a one-off `afkf`-vs-PFSimple comparison; leave it off (`false`, the default) otherwise.

### decays
One or multiple decays get can reconstructed in one run. Each decay is defined by several sub-blocks:

#### mother
- `pdg_code`
- `name` (chosen by user)
- `mass` (Unit: GeV/c^2) and `mass_sigma` (optional): if set to "-1" values from KFParticleDatabase are taken (if available) 
- `cuts`: subblock for the definition of cuts to be applied to the mother (See `at_interface/config/config_template.json` for the available JSON keys and  `Constants.hpp` for more details.)
- `save_options`: options for the saving of the reconstructed mother:

| Flag                 | Description                                                                                                                                                                                                    |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mass_constraint`    | nonlinear mass constraint is applied on the mother particle which modifies the mother's state vector and covariance matrix to assure the 4-momentum is consitent with the pdg mass                             |
| `transport_to_pv`    | mother is transported to the PV before saving                                                                                                                                                                  |
| `do_not_save_mother` | writing of the mother into the output-file (e.g. the lower-level mother in cascade decays) is surpressed. If this option is selected, it is not possible to peform the MC-matching for the upper-level mother. |

#### daughters
The 2 (3) daughters of the mother are be defined in 2 (3) subblocks (for 2-/3-body decays).
- `pdg_code` (in JSON list format):
    first entry: daugher PDG;
    second and following entries (optional): additional pdgs or background can be considered as possible daughter candidates, for background set:
    - for background: set to "1" / "-1" for pos / neg daughters
    - for tracks without TOF id: set to "2" / "-2" for pos / neg daughters (in pid-mode 2, 3, 4)
    In pid-mode 0 optional daughers "1" / "-1" for pos / neg pdgs must be added
- `cuts`:  subblock for the definition of cuts to be applied to the daughter (See `at_interface/config/config_template.json` for the available JSON keys and  `Constants.hpp` for more details.)

#### secondary_mother_cuts (3-body-decays)
Cuts on the secondary mothers (SM) of combinations of 2 daughters can have up to 3 entries (JSON list format): cuts for daughter combinations  [1 & 2, 1 & 3, 2 & 3].
E.g. for H3L->p + pi + d, a cut on the SM of p-pi can be applied.
The cuts that test the mother against the PV, e.g. chi2topo, are inverted for the SM to exclude SMs that come from the PV.
(See `at_interface/config/config_template.json` for the available JSON keys and  `Constants.hpp` for more details.)

### output_cuts
The optional `output_cuts` block contains AnalysisTree field cuts to be applied to the output candidates before storing them to disk. One can either select a certain range for a cut variable with `from` and `to` or select a specific value with `equal`. It is possible to chain several cuts. 

### event_mixing
The optional `event_mixing` block estimates the combinatorial background by reconstructing candidates from tracks that belong to different events (event mixing). It follows the track-level mixing approach of the `afkf` package: the tracks of a chosen species (`mix_pdg`) are taken from the preceding events and combined with the remaining tracks of the current event. The resulting background candidates are written into a separate `MixedCandidates` branch (cloned from the `Candidates` branch), so that they can be analysed with the same code as the signal candidates.

| Key                | Description                                                                                                                                                              |
| -------------------| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `enabled`          | Enable/disable event mixing (default: `true` when the block is present)                                                                                                  |
| `mix_pdg`          | PDG of the track species that is replaced by tracks of previous events (default: `-211`, i.e. the pion of the Lambda decay)                                             |
| `mix_buffer_size`  | Number of previous events kept in the mixing buffer (default: `5`)                                                                                                       |
| `branchname`       | Optional name of the output branch that stores the mixed-event background candidates (default: `MixedCandidates`)                                                        |
| `centrality_source`| Source of the centrality value used to bin the mixing buffer (default: `"field"`). See the centrality table below.                                                        |
| `centrality_field_name` | Name of the event-header field holding the centrality value, used when `centrality_source` is `"field"` (default: `"centrality_tracks"`). Useful for AT preprocessing variants that name the field differently. |

For each event, the `mix_pdg` tracks of the current event are discarded and replaced by the `mix_pdg` tracks of the `mix_buffer_size` preceding events, translated to the current primary vertex. The combined track set is re-id'ed to a contiguous range and reconstructed with the same decays as the signal candidates; all candidates built in this pass are combinatorial background. If MC info is available, the `generation` field of the mixed candidates is set to `-1000` to distinguish them from signal candidates.

Mixing is single-species and sign-sensitive: only daughters whose `pdg_code` exactly matches `mix_pdg` are swapped (as in `afkf`). When the config has several decays, a startup summary reports, per mother, whether it is `ALL species mixed`, `PARTIAL` (only some daughters match, so the mixed sample is a partially-mixed background), or `NOT mixed` (no daughter matches `mix_pdg`, so no mixed background is produced for it). The current event is always excluded from its own buffer (read-then-fill), so its real daughters never combine with each other; however, the mixed daughters are drawn from a flat pool of the `mix_buffer_size` preceding events in the same centrality bin, so same-species daughters can occasionally both come from one buffered event (a small cross-event leak, identical to `afkf`, but not a current-event signal).

The mixing buffer is binned by centrality (9 bins, mirroring `afkf`): an event only mixes with the last `mix_buffer_size` events of the *same* bin. `centrality_source` selects how the bin is obtained:

| `centrality_source`  | Meaning                                                                                                                                                                                           | Matches `afkf` when                                        |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------|
| `"field"` (default)  | Read the `centrality_tracks` field from the event header and bin with the `{80,70,60,50,40,30,20,10,5}` table (`CentralityAT`); bin 8 = most central (lowest centrality value), bin 0 = most peripheral. | `at_tree == "aTree"`                                       |
| `"ntracks"`          | Use the event's good-track count, binned with the `{6,12,22,38,62,95,145,212,257}` table (`CentralityNTacks`); bin 8 = most central (highest track count), bin 0 = most peripheral.               | informational only — **not** what `afkf` does for non-`aTree` |
| `"fixed8"`           | Always use bin 8 (most central); all events mix within a single buffer.                                                                                                                          | `at_tree != "aTree"` (e.g. `rTree`, `aTree2`) — **use this for the cross-check** |
| `"none"`             | No centrality binning — all events share a single buffer (bin 0).                                                                                                                                 | —                                                          |

`afkf` picks the bin depending on the input tree name. For `at_tree == "aTree"` it reads `centrality_tracks` and bins with `Centrality_AT` (`{80,70,60,50,40,30,20,10,5}`). For **non-**`aTree` inputs it **hard-codes `fCentrality = 8`** at the start of every event, so all events mix within a single bin-8 buffer regardless of their actual centrality. PFSimple mirrors this with `centrality_source: "fixed8"` for the `rTree`/`aTree2` cross-check. Using `"ntracks"` (bins by actual track count) does **not** match `afkf` for non-`aTree` inputs and leads to a ~20–30% background deficit, so it is not recommended for that comparison.

### Cascade decays
For the reconstruction of a cascade decay all consecutive decays are added into one config-file. The order of the decays need to be from last generation to
first generation.
Cascade decays with multiple stages as well as combinations of 2- and 3-body-decays can get reconstructed.

## Example Configs
Several example parameter files can be found in the folder `at_interface/configs`.

| File                                  | Particle(s)        | `pid_mode` | Comment                                                                                                                        |
|-------------------------------------- |--------------------|------------|--------------------------------------- |
| `config_template.json`                | -                  | -          | all available JSON fields listed       |
| `config_lambda_pidmode1_v1.json`      | Lambda             | `1`        | KFParticle cuts                        |
| `config_lambda_pidmode1_v2.json`      | Lambda             | `1`        | optimized cuts                         |
| `config_lambda_pidmode0.json`         | Lambda             | `0`        |                                        |
| `config_lambda_pidmode2.json`         | Lambda             | `2`        |                                        |
| `config_lambda_pidmode1_v3.json`      | Lambda             | `1`        | for machine learning applications      |
| `config_H3L_pidmode1.json`            | H3L                | `1`        | 3-body decay                           |
| `config_H3L_pidmode4.json`            | H3L                | `4`        | 3-body decay, minimum purities applied |
| `config_He5L_pidmode1.json`           | He5L               | `1`        | 3-body decay                           |
| `config_xi_pidmode1.json`             | Xi, Lambda         | `1`        | cascade decay                          |
| `config_omega_pidmode1.json`          | Omega, Lambda      | `1`        | cascade decay                          |
| `config_omegastar_pidmode1.json`      | Omega*, Xi, Lambda | `1`        | cascade decay, 3-body decay            |
| `config_lambda_k0short_pidmode1.json` | Lambda, K0Short    | `1`        |                                        |
| `config_lambda_mix_afkf_xcheck.json`  | Lambda             | `1`        | event mixing, reproduces the `afkf` config `atkf_runtime_config_mix_pion.json` for cross-validation (incl. `mc_match_weight_min: 0.7` to mirror afkf's hardcoded track selection, `accept_fragments: true` and `chi2_prim_afkf: true` to mirror afkf's fragment acceptance and chi2-to-PV) |



## First run
Each time before running the prepared executable you should set the environment variables to let your system know where to find libraries:
```
source /path-to-root/install/bin/thisroot.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/external/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-analysistree-installation/lib/ #can be ommitted if PFSimple_BUNDLED_AT was set to ON
```

The reconstruction will be executed with the following command:
```
at_interface/main_json filelist.txt config.json [--output <file>] [--plain-output <file>] [--set <key>=<value>]
```
The binary accepts 2 positional arguments: The path of the input filelist.txt file and the config file in JSON format. It is also possible to modify the output file paths with the `--output` and `--plain-output` optional arguments as well as to override keys in the JSON config file by adding keys like `--set decays[0].mother.cuts.LdL=4.0`.

### filelist.txt
filelist.txt must be a text file with names (including paths) to AnalysisTree
files which you want to analyze with PFSimple. Each file name should
be on the next line, and the last symbol should also be a switch to
next line. The example of filelist is in the source directory of
PFSimple (filelist_example.txt)
