use std::io;
use tracing::Level;
use tracing_subscriber::{filter, fmt, prelude::*, EnvFilter};

pub fn start_skim_tracing_subscriber() -> () {
    // Create a layer that logs to stdout
    let stdout_log = fmt::layer();

    // Get the stdout logging filter level from the RUST_LOG environment variable
    //   - INFO messages are always logged to stdout
    //   - If RUST_LOG=debug, DEBUG messages will also be logged to stdout
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));

    // Create a layer that logs to stderr
    let stderr_log = fmt::layer().with_writer(io::stderr);

    // Combine the layers into a registry (a subscriber)
    // For the stdout log, make sure INFO is always logged and DEBUG is also logged if requested
    // For the stderr log, log WARN and ERROR
    let tracer_registry = tracing_subscriber::registry()
        .with(
            stdout_log
                .with_filter(env_filter)
                .with_filter(filter::filter_fn(|metadata| {
                    *metadata.level() >= Level::INFO
                })),
        )
        .with(stderr_log.with_filter(filter::LevelFilter::WARN));

    // Initialize the subscriber
    tracer_registry.init()
}
