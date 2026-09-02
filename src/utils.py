"""
Utility functions for single-cell preprocessing pipeline.
Handles logging, configuration, command execution, and reporting.
"""

import os
import sys
import subprocess
import shutil
import yaml
from pathlib import Path
from typing import Dict, Any, List, Optional
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

def load_config(config_path: str = "config/pipeline_config.yaml") -> Dict[str, Any]:
    """Load and validate pipeline configuration YAML file."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with open(path, "r") as f:
        config = yaml.safe_load(f)
    return config

def ensure_dir(dir_path: str | Path) -> Path:
    """Ensure directory exists and return Path object."""
    p = Path(dir_path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def run_cmd(
    cmd: List[str] | str,
    desc: Optional[str] = None,
    check: bool = True,
    capture_output: bool = True
) -> subprocess.CompletedProcess:
    """
    Run an external command with clean logging and error reporting.
    """
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
        shell = False
    else:
        cmd_str = cmd
        shell = True
        
    if desc:
        console.print(f"[bold cyan]▶ {desc}[/bold cyan]: [dim]{cmd_str}[/dim]")
        
    try:
        res = subprocess.run(
            cmd if not shell else cmd_str,
            shell=shell,
            check=check,
            capture_output=capture_output,
            text=True
        )
        return res
    except subprocess.CalledProcessError as e:
        console.print(f"[bold red]✘ Command failed ({e.returncode}):[/bold red] {cmd_str}")
        if e.stderr:
            console.print(f"[red]Error Output:[/red]\n{e.stderr.strip()[:1000]}")
        raise e

def print_summary_table(title: str, metrics: Dict[str, Any]):
    """Display a beautifully formatted summary table in console."""
    table = Table(title=title, header_style="bold magenta", border_style="cyan")
    table.add_column("Metric", style="bold white")
    table.add_column("Value", style="green", justify="right")
    
    for k, v in metrics.items():
        if isinstance(v, float):
            val_str = f"{v:.4f}" if abs(v) < 1 else f"{v:,.2f}"
            if "fraction" in k.lower() or "percent" in k.lower() or "rate" in k.lower() or "q30" in k.lower():
                val_str = f"{v*100:.2f}%" if v <= 1.0 else f"{v:.2f}%"
        elif isinstance(v, int):
            val_str = f"{v:,}"
        else:
            val_str = str(v)
        table.add_row(k.replace("_", " ").title(), val_str)
        
    console.print(table)

def execute_notebook(notebook_path: str | Path, output_path: Optional[str | Path] = None) -> bool:
    """Execute a Jupyter notebook programmatically to verify correctness."""
    import nbformat
    from nbclient import NotebookClient
    
    nb_path = Path(notebook_path)
    if not nb_path.exists():
        console.print(f"[red]Notebook not found: {nb_path}[/red]")
        return False
        
    console.print(f"[bold yellow]Executing notebook:[/bold yellow] {nb_path.name}")
    with open(nb_path, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)
        
    client = NotebookClient(nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(nb_path.parent.resolve())}})
    try:
        client.execute()
        if output_path:
            out_p = Path(output_path)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            with open(out_p, "w", encoding="utf-8") as f:
                nbformat.write(nb, f)
        console.print(f"[bold green]✔ Notebook executed successfully:[/bold green] {nb_path.name}")
        return True
    except Exception as e:
        console.print(f"[bold red]✘ Error executing notebook {nb_path.name}:[/bold red] {e}")
        return False

if __name__ == "__main__":
    if "--test-notebooks" in sys.argv:
        notebooks = sorted(Path("notebooks").glob("*.ipynb"))
        all_passed = True
        for nb in notebooks:
            if not execute_notebook(nb):
                all_passed = False
        sys.exit(0 if all_passed else 1)
