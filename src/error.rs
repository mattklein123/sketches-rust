use std::io;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Invalid argument: {0}")]
    InvalidArgument(&'static str),
    #[error("IO error: {0}")]
    IoError(io::ErrorKind),
}
