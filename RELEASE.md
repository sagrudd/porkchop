# Release Checklist

Target version: 0.2.29

## Pre-flight
- [ ] Update `CHANGELOG.md` with notable changes.
- [ ] Ensure `Cargo.toml` and README **Version** agree.
- [ ] `cargo fmt --all` (clean format).
- [ ] `cargo clippy --all-targets --all-features -D warnings` (no warnings).
- [ ] `cargo test --all-features` (all green).
- [ ] Verify MSRV job is green in CI (1.82.0).

## Tag & Publish
- [ ] Create tag `v0.2.29`.
- [ ] `cargo publish` (or release binaries if applicable).
- [ ] Create GitHub release with changelog excerpt.

_Generated on 2025-10-23._
