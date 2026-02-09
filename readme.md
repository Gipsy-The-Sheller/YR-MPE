# YR-MPE - Molecular Phylogenetics & Evolution Plugins for YRTools

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-GPL-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/Gipsy-The-Sheller/YR-MPE)

<img src="/YR_MPE/icons/yr-mpea.svg" width="auto" height="50"/>

> [!NOTE]
> 
> If you find this repository useful for your research, please consider giving it a star! It helps me gauge interest. For any questions, please open an issue.
>
> If you want a specific function, please describe it in detail in a new issue. I will consider implementing it in priority.
>
> YR-MPE temporarily have no publications. If you need to cite it, please cite the repository.

## What is YR-MPE?

`YR-MPE` is originally a collection of plugins for `YRTools`, which is a plugin-based bioinformatic platform. `YR-MPE` is for its molecular phylogenetics and evolution analysis functions.

While other software also provides good experience, an inherent conflict is the gap between these integrative software's good experience and most up-to-date functions, since computational molecular evolution is developed at a rapid pace.

Owing to a good architecture, the source codes of `YR-MPE` do not need any compilation or packaging, even without a global python environment. Thus, every new commit can be a new release, supporting speedy updates directly from Github repository. Also, it is open-sourced, which means you can modify it as you wish, and community can develop new functions in collaboration. As a result, I think `YR-MPE` is capable to keep up with the development of the discipline.

## How to use it?

You can invoke all or part of `YR-MPE` from any UI software (or even a python script) in `PyQt5` ecology.

The most convenient way is to use it in `YRTools` [Github](https://github.com/Gipsy-The-Sheller/YRTools) [Gitee-coming out soon](). You can get the up-to-date version direct from repository on click with its package manager `YR-Pacman`.

The wiki of `YR-MPE` provides detailed introduction and usage. You can access its [Github](https://github.com/Gipsy-The-Sheller/YR-MPE/wiki) or [Gitee](https://gitee.com/ZJXMolls/YR-MPE/wikis) source online.

## Changelogs

2026.1.9 - v0.2.0 - Bayesian inference plugins (MrBayes-MPI-BEAGLE3, PhyloBayes-MPI), fixed IcyTree bugs, added console plugin.