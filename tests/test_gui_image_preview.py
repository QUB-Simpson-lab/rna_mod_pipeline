from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication, QTreeWidgetItem

from rna_mod_pipeline.gui.image_preview import ImagePreview
from rna_mod_pipeline.gui.results_page import ResultsPage


class ImagePreviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.widgets = []

    def tearDown(self) -> None:
        for widget in self.widgets:
            widget.close()
            widget.deleteLater()
        self.app.processEvents()
        self.temporary.cleanup()

    def make_image(self, width=2382, height=821) -> Path:
        path = self.root / f"plot {width}x{height}.png"
        image = QImage(width, height, QImage.Format.Format_RGB32)
        image.fill(Qt.GlobalColor.white)
        self.assertTrue(image.save(str(path)))
        return path

    def viewer(self) -> ImagePreview:
        viewer = ImagePreview()
        self.widgets.append(viewer)
        viewer.resize(700, 500)
        viewer.show()
        self.app.processEvents()
        return viewer

    def assert_fits(self, viewer: ImagePreview) -> None:
        self.app.processEvents()
        mapped = viewer.view.mapFromScene(viewer.scene.sceneRect()).boundingRect()
        viewport = viewer.view.viewport().rect()
        self.assertTrue(viewport.contains(mapped), (viewport, mapped))
        transform = viewer.view.transform()
        self.assertAlmostEqual(transform.m11(), transform.m22())
        self.assertTrue(viewer.fit_mode)

    def test_wide_and_tall_images_fit_without_cropping(self) -> None:
        viewer = self.viewer()
        for width, height in ((2382, 821), (2198, 7955), (6004, 1639)):
            with self.subTest(dimensions=(width, height)):
                path = self.make_image(width, height)
                before = path.read_bytes()
                viewer.set_image(path)
                self.assertEqual(viewer.path, path)
                self.assert_fits(viewer)
                self.assertEqual(path.read_bytes(), before)

    def test_fit_tracks_window_resize(self) -> None:
        viewer = self.viewer()
        viewer.set_image(self.make_image())
        old_scale = viewer.scale_factor
        viewer.resize(450, 350)
        self.assert_fits(viewer)
        self.assertLess(viewer.scale_factor, old_scale)
        viewer.resize(950, 750)
        self.assert_fits(viewer)
        self.assertGreater(viewer.scale_factor, old_scale)

    def test_actual_size_scroll_zoom_and_return_to_fit(self) -> None:
        viewer = self.viewer()
        viewer.set_image(self.make_image())
        viewer.actual_button.click()
        self.app.processEvents()
        self.assertEqual(viewer.scale_factor, 1.0)
        self.assertGreater(viewer.view.horizontalScrollBar().maximum(), 0)
        self.assertGreater(viewer.view.verticalScrollBar().maximum(), 0)
        viewer.plus_button.click()
        self.assertEqual(viewer.scale_factor, 1.25)
        viewer.resize(800, 600)
        self.app.processEvents()
        self.assertEqual(viewer.scale_factor, 1.25)
        viewer.minus_button.click()
        self.assertEqual(viewer.scale_factor, 1.0)
        viewer.fit_button.click()
        self.assert_fits(viewer)

    def test_corrupt_file_clears_old_image_and_disables_controls(self) -> None:
        viewer = self.viewer()
        viewer.set_image(self.make_image())
        broken = self.root / "broken.png"
        broken.write_bytes(b"not an image")
        viewer.set_image(broken)
        self.assertIsNone(viewer.path)
        self.assertFalse(viewer.scene.items())
        self.assertIn("Cannot preview broken.png", viewer.caption.text())
        self.assertTrue(all(not button.isEnabled() for button in viewer.controls))

    def test_open_original_uses_local_file_url(self) -> None:
        viewer = self.viewer()
        image = self.make_image()
        viewer.set_image(image)
        with patch("rna_mod_pipeline.gui.image_preview.QDesktopServices.openUrl") as opened:
            viewer.open_button.click()
        self.assertEqual(opened.call_args.args[0].toLocalFile(), str(image))

    def page(self) -> ResultsPage:
        page = ResultsPage(self.root)
        self.widgets.append(page)
        page.resize(1200, 800)
        page.show()
        self.app.processEvents()
        return page

    def choose(self, page: ResultsPage, path: Path) -> None:
        item = QTreeWidgetItem([path.name])
        item.setData(0, Qt.ItemDataRole.UserRole, str(path))
        page.tree.addTopLevelItem(item)
        page.tree.setCurrentItem(item)
        self.app.processEvents()

    def test_images_use_full_panel_and_text_switches_cleanly(self) -> None:
        page = self.page()
        self.choose(page, self.make_image())
        self.assertFalse(page.text_preview.isVisible())
        self.assertTrue(page.image.isVisible())
        self.assertGreater(page.image.view.viewport().height(), 500)
        self.assert_fits(page.image)
        table = self.root / "example.tsv"
        table.write_text("RBP\tOR\nDEMO\t9\n")
        self.choose(page, table)
        self.assertTrue(page.text_preview.isVisible())
        self.assertFalse(page.image.isVisible())
        self.assertIsNone(page.image.path)
        self.assertIn("DEMO", page.text_preview.toPlainText())
        self.choose(page, self.make_image(500, 1200))
        self.assertFalse(page.text_preview.isVisible())
        self.assert_fits(page.image)

    def test_missing_image_shows_error_instead_of_previous_plot(self) -> None:
        page = self.page()
        self.choose(page, self.make_image())
        self.choose(page, self.root / "missing.png")
        self.assertIsNone(page.image.path)
        self.assertFalse(page.image.isVisible())
        self.assertIn("missing", page.text_preview.toPlainText())

    def test_new_manifest_clears_previous_preview(self) -> None:
        page = self.page()
        self.choose(page, self.make_image())
        manifest = self.root / "run_manifest.json"
        manifest.write_text(json.dumps({"workflow": "new", "outputs": []}))
        page.load_output(manifest)
        self.assertEqual(page.view.workflow, "new")
        self.assertIsNone(page.image.path)
        self.assertEqual(page.tree.topLevelItemCount(), 0)

    def test_invalid_manifest_does_not_leave_stale_outputs(self) -> None:
        page = self.page()
        self.choose(page, self.make_image())
        with patch("rna_mod_pipeline.gui.results_page.QMessageBox.critical") as error:
            page.load_output(self.root / "missing.json")
        error.assert_called_once()
        self.assertIsNone(page.view)
        self.assertIsNone(page.image.path)
        self.assertEqual(page.tree.topLevelItemCount(), 0)


if __name__ == "__main__":
    unittest.main()
