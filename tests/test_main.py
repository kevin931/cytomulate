from cytomulate.__main__ import main, _Cli
import cytomulate

import argparse
from io import StringIO
import sys
import pytest
from typing import List


class TestCli():
    
    @classmethod
    def setup_class(cls):
        cls.arguments_parser: _Cli = _Cli()
        
        
    def test_parser_type(self):
        assert isinstance(self.arguments_parser.parser, argparse.ArgumentParser)
        
    
    # TODO: test_parse namespace
    
    
    @pytest.mark.parametrize("arguments",
                            [["-h"], ["--version"]])
    def test_version_system_exit(self, arguments):
        try:
            self.arguments_parser.parse(arguments)
        except SystemExit:
            assert True
        else:
            assert False
       
            
    @pytest.mark.parametrize("arguments,expected",
                            [(["-h"], "show this help message"),
                            (["--version"], cytomulate.__version__)]
                            )
    def test_cli_message(self, arguments: List[str], expected: str):
        screen_stdout = sys.stdout
        string_stdout = StringIO()
        sys.stdout = string_stdout
        
        try: 
            self.arguments_parser.parse(arguments)
        except SystemExit:
            output = string_stdout.getvalue()
            sys.stdout = screen_stdout
            assert expected in output
        else:
            sys.stdout = screen_stdout
            assert False
            