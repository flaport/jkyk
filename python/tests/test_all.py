import pytest
import jkyk


def test_sum_as_string():
    assert jkyk.sum_as_string(1, 1) == "2"
