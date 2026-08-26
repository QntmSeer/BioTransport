//! CLI interface for the biotransport high-performance solver engine.

use clap::{Parser, Subcommand};
use std::fs;
use std::path::PathBuf;

use biotransport::{
    simulate_tff_filtration, MdBridgeParameters, ProcessOperatingConditions,
};

#[derive(Parser)]
#[command(name = "biotransport-cli")]
#[command(about = "Continuum membrane transport & fouling solver parameterized by CG-MD data", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run a TFF filtration simulation using bridged MD parameters
    Run {
        /// Path to input JSON with MD parameters
        #[arg(short, long)]
        params: PathBuf,

        /// Transmembrane pressure in Pascals (default: 150000 Pa = 1.5 bar)
        #[arg(long, default_value_t = 150_000.0)]
        tmp: f64,

        /// Bulk feed concentration in g/L (default: 10.0 g/L)
        #[arg(long, default_value_t = 10.0)]
        bulk_conc: f64,

        /// Total filtration time in seconds (default: 3600 s)
        #[arg(short, long, default_value_t = 3600.0)]
        time: f64,

        /// Path to write output results JSON
        #[arg(short, long)]
        out: Option<PathBuf>,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Run {
            params,
            tmp,
            bulk_conc,
            time,
            out,
        } => {
            println!("Reading MD parameters from: {:?}", params);
            let content = fs::read_to_string(&params)?;
            let md_params: MdBridgeParameters = serde_json::from_str(&content)?;

            let ops = ProcessOperatingConditions {
                transmembrane_pressure_pa: tmp,
                bulk_concentration_g_l: bulk_conc,
                total_time_s: time,
                ..Default::default()
            };

            println!("Running continuum filtration solver...");
            let summary = simulate_tff_filtration(&md_params, &ops);

            println!("================ Simulation Results ================");
            println!("Initial Flux: {:.2} LMH", summary.initial_flux_lmh);
            println!("Final Flux:   {:.2} LMH", summary.final_flux_lmh);
            println!("Flux Decline: {:.1}%", summary.flux_decline_percentage);
            println!("Total Permeate: {:.2} L/m^2", summary.total_permeate_collected_l_m2);
            println!("Max Wall Conc:  {:.1} g/L", summary.maximum_wall_concentration_g_l);
            println!("Specific Cake Resistance: {:.2e} m/kg", summary.specific_cake_resistance_m_kg);
            println!("====================================================");

            if let Some(out_path) = out {
                let serialized = serde_json::to_string_pretty(&summary)?;
                fs::write(&out_path, serialized)?;
                println!("Saved full trajectory to {:?}", out_path);
            }
        }
    }

    Ok(())
}
