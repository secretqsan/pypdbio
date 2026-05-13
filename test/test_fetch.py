import os
import tempfile
import unittest
from pypdbio import fetch

class TestFetch(unittest.TestCase):
    def test_fetch_can_download_sample_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "1aki.pdb")
            fetch("1aki", path=out_path)
            self.assertTrue(os.path.exists(out_path))

if __name__ == "__main__":
    unittest.main()