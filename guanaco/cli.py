#!/usr/bin/env python3
"""Command-line interface for GUANACO visualization tool."""

import argparse
import os
from pathlib import Path

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="GUANACO: Interactive visualization tool for single-cell and genome browser data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-c", "--config",
        type=str,
        default=Path("guanaco.json"),
        help="Name of configuration JSON file(relative to --data-dir)"
    )
    
    parser.add_argument(
        "-d", "--data-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory containing AnnData files referenced in config"
    )
    
    parser.add_argument(
        "-p", "--port",
        type=int,
        default=4399,
        help="Port to run the Dash server on"
    )
    
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Host to run the Dash server on"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run server in debug mode"
    )
    
    parser.add_argument(
        "--max-cells",
        type=int,
        default=8000,
        help="Maximum number of cells to load per dataset"
    )
    
    
    args = parser.parse_args()
    
    # Validate config file exists
    config_path = args.data_dir / args.config
    if not config_path.exists():
        parser.error(f"Config file '{config_path}' not found.")

    
    # Validate data directory exists
    if not args.data_dir.exists():
        parser.error(f"Data directory '{args.data_dir}' not found.")   
    return args


def main():
    """Main entry point for GUANACO CLI."""
    args = parse_args()
    
    # Set environment variables from CLI args
    os.environ['GUANACO_DATA_DIR'] = str(args.data_dir)
    config_path = (args.data_dir / args.config).resolve()
    os.environ['GUANACO_CONFIG'] = str(config_path)
    os.environ['GUANACO_MAX_CELLS'] = str(args.max_cells)
    
    print("Starting GUANACO server...")
    print(f"Config file: {args.config}")
    print(f"Data directory: {args.data_dir}")
    print(f"Server will run on {args.host}:{args.port}")
    print(f"Debug mode: {args.debug}")
    from guanaco.main import app
    app.run_server(host=args.host, debug=args.debug, port=args.port)

if __name__ == "__main__":
    main()
