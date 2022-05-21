from cytomulate import __main__
import cytomulate
import pytest

import sys
from io import StringIO


def test_parser_version():
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args(["--version"])
    except SystemExit:
        output = string_stdout.getvalue()
        expected = cytomulate.__version__ + "\n"
        assert output == expected
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False
        

@pytest.mark.parametrize("cmdarg", ["--help", "-h"])
def test_parser_help(cmdarg):
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args([cmdarg])
    except SystemExit:
        output = string_stdout.getvalue()
        assert "usage: " in output
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False