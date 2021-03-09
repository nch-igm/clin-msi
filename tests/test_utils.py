from unittest import TestCase

from clin_msi import utils


class TestRepeatFinder(TestCase):
    
    def test_repeat_finder(self):
        assert list(utils.repeat_finder("blablabla")) == [('bla', 3)]
        assert list(utils.repeat_finder("rablabla")) == [('abl', 2)]
        assert list(utils.repeat_finder("aaaaa")) == [('a', 5)]
        assert list(utils.repeat_finder("aaaaablablabla")) == [('a', 5), ('bla', 3)]
