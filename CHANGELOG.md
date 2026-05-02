# Changelog

## [1.4.0](https://github.com/PathoGenOmics-Lab/VIPERA/compare/v1.3.1...v1.4.0) (2026-05-02)


### Features

* add low-confidence site filtering ([f447478](https://github.com/PathoGenOmics-Lab/VIPERA/commit/f447478fe409e01e64c7a9d077a18c8518f460aa))
* remove demix module from rule core ([#30](https://github.com/PathoGenOmics-Lab/VIPERA/issues/30)) ([cb01079](https://github.com/PathoGenOmics-Lab/VIPERA/commit/cb0107939b82a2d7cc0dc94ca69448f9e6cf5b83))
* use site quality filter to calculate evolutionary rate ([#33](https://github.com/PathoGenOmics-Lab/VIPERA/issues/33)) ([29d75e5](https://github.com/PathoGenOmics-Lab/VIPERA/commit/29d75e5e30d93b7c1048d93cda99dabc5a21690f))


### Bug Fixes

* apply negative branch length correction iteratively ([#45](https://github.com/PathoGenOmics-Lab/VIPERA/issues/45)) ([76b2957](https://github.com/PathoGenOmics-Lab/VIPERA/commit/76b29578710a9f75a443f2665605eda73a12519b))
* catch errors on AF-time correlation test ([#39](https://github.com/PathoGenOmics-Lab/VIPERA/issues/39)) ([8025462](https://github.com/PathoGenOmics-Lab/VIPERA/commit/8025462be3d0710f34216307a739e79d00dd16ac))
* collapse multi-site variant names ([#42](https://github.com/PathoGenOmics-Lab/VIPERA/issues/42)) ([b9fee30](https://github.com/PathoGenOmics-Lab/VIPERA/commit/b9fee305132448cd45532b68065667cfcd112e1f))
* consider masked sites for dN and dS calculation ([#38](https://github.com/PathoGenOmics-Lab/VIPERA/issues/38)) ([ce47836](https://github.com/PathoGenOmics-Lab/VIPERA/commit/ce478362f6a34997add99932deaf1403b731d769))
* deduplicate variants and handle missing annotation ([#32](https://github.com/PathoGenOmics-Lab/VIPERA/issues/32)) ([c70dfe0](https://github.com/PathoGenOmics-Lab/VIPERA/commit/c70dfe0a606705807c3f912a21a89338a9f88625))
* handle empty extracted VCF fields ([#43](https://github.com/PathoGenOmics-Lab/VIPERA/issues/43)) ([a90ca5b](https://github.com/PathoGenOmics-Lab/VIPERA/commit/a90ca5bb5859e907aa3ce1c7630de93982aa8637))
* init missing columns ([#37](https://github.com/PathoGenOmics-Lab/VIPERA/issues/37)) ([fa5f624](https://github.com/PathoGenOmics-Lab/VIPERA/commit/fa5f6240f4eac915eda1a61c53fbce0db8370732))
* pin setuptools to avoid deprecation error in Snakemake 7 ([395696c](https://github.com/PathoGenOmics-Lab/VIPERA/commit/395696c8f1231d8e94c5eaf330a6dda59045d210))
* pin setuptools version to avoid deprecation error ([#31](https://github.com/PathoGenOmics-Lab/VIPERA/issues/31)) ([395696c](https://github.com/PathoGenOmics-Lab/VIPERA/commit/395696c8f1231d8e94c5eaf330a6dda59045d210))

## [1.3.1](https://github.com/PathoGenOmics-Lab/VIPERA/compare/v1.3.0...v1.3.1) (2025-11-24)


### Bug Fixes

* force reading VCF fields as character ([18e5004](https://github.com/PathoGenOmics-Lab/VIPERA/commit/18e5004b9a3d8bfe0975109d2832a921d047930c))
