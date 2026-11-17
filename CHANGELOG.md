# Changelog

## 0.1.3 - 2026-11-17

### Fixed

- sse: fix opstring estimators failing on Julia nightly, which passes `Core.TypeEgal` instead of `Type` to generated functions.
- make the deterministic S=1/2 case actually trigger (thanks @SufficientlySmooth!)

## 0.1.2 - 2025-12-04

### Changed

- update Carlo compat 

## 0.1.1 - 2024-11-21

### Added
- `ClusterModel` for simulating other models in clustered bases (e.g. fully-frustrated bilayer)
### Changed
- `MagnetizationEstimator`: add `Tag` parameter to allow different kinds of magnetizations for the same model
### Fixed
- magnet: fix typo in the Hamiltonian that changed the meaning of `d`.
- some performance optimizations
