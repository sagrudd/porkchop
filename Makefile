    .PHONY: fmt clippy check test all

    fmt:
	cargo fmt --all

    clippy:
	cargo clippy --all-targets --all-features -D warnings

    check:
	cargo check --all-targets --all-features

    test:
	cargo test --all-features

    all: fmt clippy check test
