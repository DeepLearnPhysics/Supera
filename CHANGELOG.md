# Changelog

All notable changes to Supera will be documented in this file.

## [1.0.0] - 2026-06-28

### Added
- Began explicit Supera version tracking.

### Changed
- Optimized `SuperaSpacePoint` voxel accumulation and skipped filling configured dropped outputs.
- Adapted `SuperaOptical` to the LArCV2 `Flash` volume-id field.
- Preserved legacy SBND production FHiCL settings.

### Fixed
- Fixed the `SuperaCRT` include guard so it no longer collides with `SuperaOptical`.
