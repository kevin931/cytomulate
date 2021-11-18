from cytomulate import __main__
import cytomulate

import argparse
from io import StringIO
import sys
import pytest
from typing import List


class TestCli():    
        
    def test_parser_type(self):
        assert isinstance(__main__.parser, argparse.ArgumentParser)
        
    
    # TODO: test_parse namespace
    
    
    @pytest.mark.parametrize("arguments",
                            [["-h"], ["--version"]])
    def test_version_system_exit(self, arguments):
        try:
            __main__.parser.parse_args(arguments)
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
            __main__.parser.parse_args(arguments)
        except SystemExit:
            output = string_stdout.getvalue()
            sys.stdout = screen_stdout
            assert expected in output
        else:
            sys.stdout = screen_stdout
            assert False
            