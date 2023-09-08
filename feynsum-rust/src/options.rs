use std::path::PathBuf;
use structopt::StructOpt;

use crate::gate_scheduler::GateSchedulingPolicy;

#[derive(Debug, StructOpt)]
#[structopt(name = "feynsum", about = "Feynsum quantum simulator")]
pub struct Options {
    #[structopt(
        parse(from_os_str),
        name = "input",
        short = "i",
        long = "input",
        help = "path to the input qasm file"
    )]
    pub input: PathBuf,

    #[structopt(
        parse(from_os_str),
        name = "output",
        short = "o",
        long = "output",
        help = "output file path to dump densities to. if not specified, densities are not printed"
    )]
    pub output: Option<PathBuf>,

    #[structopt(
        name = "gate scheduling policy",
        long = "scheduler",
        short = "s",
        default_value = "naive",
        help = "gate scheduling policy to use"
    )]
    pub gate_schduling_policy: GateSchedulingPolicy,

    #[structopt(
        name = "dense threshold",
        long = "dense-threshold",
        short = "dt",
        default_value = "0.25"
    )]
    pub dense_threshold: f64,

    #[structopt(
        name = "pull threshold",
        long = "pull-threshold",
        short = "pt",
        default_value = "0.8"
    )]
    pub pull_threshold: f64,
}
