"""Unit tests for pypdbio.utils chain-id helpers."""
# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
import unittest

from pypdbio.utils import chain_id_of_index, index_of_chain_id, next_chain_id


class TestUtils(unittest.TestCase):
    def test_next_chain_id(self):
        self.assertEqual(next_chain_id("A", []), "B")
        self.assertEqual(next_chain_id("A", ["B"]), "C")
        self.assertEqual(next_chain_id("A", ["B", "C", "D"]), "E")

    def test_index_of_chain_id(self):
        self.assertEqual(index_of_chain_id("A", []), 0)
        self.assertEqual(index_of_chain_id("Z", ['A', 'B', 'C']), 22)

    def test_chain_id_of_index(self):
        self.assertEqual(chain_id_of_index(0, []), "A")
        self.assertEqual(chain_id_of_index(1, ['A', 'B', 'C']), "D")

if __name__ == "__main__":
    unittest.main()
