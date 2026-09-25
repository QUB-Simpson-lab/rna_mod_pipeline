from __future__ import annotations

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from rna_mod_pipeline.gui.adapters import FieldSpec
from rna_mod_pipeline.gui.form_controls import DynamicForm


class DynamicFormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])

    def test_advanced_fields_are_collapsed_but_keep_values(self) -> None:
        form = DynamicForm()
        form.populate(
            (
                FieldSpec("basic", "Basic", default="shown"),
                FieldSpec("advanced", "Advanced", default="retained", advanced=True),
            ),
            lambda *_args: None,
        )
        self.assertTrue(form.advanced_button.isVisibleTo(form))
        self.assertFalse(form.advanced_widget.isVisible())
        self.assertEqual(
            form.values(),
            {"basic": "shown", "advanced": "retained"},
        )
        form.advanced_button.setChecked(True)
        self.assertTrue(form.advanced_widget.isVisibleTo(form))

    def test_empty_optional_list_is_returned_as_empty(self) -> None:
        form = DynamicForm()
        form.populate(
            (FieldSpec("optional", "Optional", kind="list", default=""),),
            lambda *_args: None,
        )
        self.assertEqual(form.values()["optional"], [])

    def test_numeric_fields_are_typed_and_invalid_text_is_rejected(self) -> None:
        form = DynamicForm()
        form.populate(
            (
                FieldSpec("count", "Count", kind="integer", default=20),
                FieldSpec("threshold", "Threshold", kind="float", default=0.05),
            ),
            lambda *_args: None,
        )
        self.assertEqual(form.values(), {"count": 20, "threshold": 0.05})
        form.fields["count"][1].setText("not-a-number")
        with self.assertRaisesRegex(ValueError, "Count must be an integer"):
            form.values()


if __name__ == "__main__":
    unittest.main()
