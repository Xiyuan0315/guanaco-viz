#!/usr/bin/env python3
"""Command-line interface for GUANACO visualization tool."""

import argparse
import sys
from pathlib import Path


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="GUANACO: Interactive visualization tool for single-cell and genome browser data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=Path("guanaco_v2.json"),
        help="Path to configuration JSON file"
    )
    
    parser.add_argument(
        "-d", "--data-dir",
        type=Path,
        default=Path("custom"),
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
    
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for cell subsampling"
    )
    
    args = parser.parse_args()
    
    # Validate config file exists
    if not args.config.exists():
        print(f"Error: Config file '{args.config}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Validate data directory exists
    if not args.data_dir.exists():
        print(f"Error: Data directory '{args.data_dir}' not found.", file=sys.stderr)
        sys.exit(1)
    
    return args


def main():
    """Main entry point for GUANACO CLI."""
    args = parse_args()
    
    # Import here to avoid circular imports
    import os
    import sys
    
    # Set environment variables from CLI args
    os.environ['GUANACO_CONFIG'] = str(args.config)
    os.environ['GUANACO_DATA_DIR'] = str(args.data_dir)
    os.environ['GUANACO_MAX_CELLS'] = str(args.max_cells)
    if args.seed is not None:
        os.environ['GUANACO_SEED'] = str(args.seed)
    
    print(f"Starting GUANACO server...")
    print(f"Config file: {args.config}")
    print(f"Data directory: {args.data_dir}")
    print(f"Server will run on {args.host}:{args.port}")
    print(f"Debug mode: {args.debug}")
    
    # Import and run the main app
    from guanaco.main import app
    app.run_server(host=args.host, debug=args.debug, port=args.port)


if __name__ == "__main__":
    main()