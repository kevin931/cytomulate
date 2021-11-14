import cytomulate
import argparse

from typing import Optional, Sequence


def main(args: argparse.Namespace):
    pass 


class _Cli():
    
    def __init__(self) -> None:
        self.parser: "argparse.ArgumentParser" = argparse.ArgumentParser(description="cytomulate: CyTOF Simulation")
        self.parser.add_argument("--version", action="version", version=cytomulate.__version__)
    
    
    def parse(self, args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        return self.parser.parse_args(args)


if __name__ == "__main__":
    cmdargs: argparse.Namespace = _Cli().parse()
    main(args = cmdargs)