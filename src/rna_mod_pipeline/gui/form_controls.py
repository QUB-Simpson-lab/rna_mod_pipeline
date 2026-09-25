from __future__ import annotations

from typing import Any, Callable, Iterable

from PySide6.QtCore import Qt
from PySide6.QtGui import QDoubleValidator, QIntValidator
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QLineEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .adapters import FieldSpec
from .widgets import PathEdit


class DynamicForm(QWidget):
    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        self.basic_widget = QWidget()
        self.basic_layout = QFormLayout(self.basic_widget)
        self.advanced_button = QToolButton()
        self.advanced_button.setText("Advanced settings")
        self.advanced_button.setCheckable(True)
        self.advanced_button.setChecked(False)
        self.advanced_button.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon
        )
        self.advanced_button.setArrowType(Qt.ArrowType.RightArrow)
        self.advanced_widget = QWidget()
        self.advanced_layout = QFormLayout(self.advanced_widget)
        self.advanced_widget.setVisible(False)
        self.advanced_button.toggled.connect(self._toggle_advanced)
        outer.addWidget(self.basic_widget)
        outer.addWidget(self.advanced_button)
        outer.addWidget(self.advanced_widget)
        self.fields: dict[str, tuple[FieldSpec, QWidget]] = {}

    def populate(
        self,
        fields: Iterable[FieldSpec],
        on_change: Callable[..., None],
    ) -> None:
        for layout in (self.basic_layout, self.advanced_layout):
            while layout.rowCount():
                layout.removeRow(0)
        self.fields.clear()
        advanced_count = 0
        for field in fields:
            widget = self._widget(field, on_change)
            label = field.label + (" *" if field.required else "")
            layout = self.advanced_layout if field.advanced else self.basic_layout
            layout.addRow(label, widget)
            advanced_count += int(field.advanced)
            widget.setToolTip(field.help_text)
            self.fields[field.key] = (field, widget)
        self.advanced_button.setVisible(bool(advanced_count))
        self.advanced_button.setChecked(False)
        self.advanced_widget.setVisible(False)

    def _toggle_advanced(self, expanded: bool) -> None:
        self.advanced_button.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self.advanced_widget.setVisible(expanded)

    def values(self) -> dict[str, Any]:
        result = {}
        for key, (field, widget) in self.fields.items():
            if isinstance(widget, QComboBox):
                value: Any = widget.currentText()
            elif isinstance(widget, QCheckBox):
                value = widget.isChecked()
            elif isinstance(widget, PathEdit):
                value = widget.text()
            else:
                value = widget.text().strip()
            kind = field.kind.lower()
            if kind in {"list", "multi"} and isinstance(value, str):
                value = [item.strip() for item in value.split(",") if item.strip()]
            elif kind in {"integer", "int"} and isinstance(value, str):
                try:
                    value = int(value)
                except ValueError as exc:
                    raise ValueError(f"{field.label} must be an integer") from exc
            elif kind in {"float", "number"} and isinstance(value, str):
                try:
                    value = float(value)
                except ValueError as exc:
                    raise ValueError(f"{field.label} must be a number") from exc
            result[key] = value
        return result

    def missing_required(self) -> list[str]:
        missing = []
        for field, widget in self.fields.values():
            if not field.required:
                continue
            if isinstance(widget, PathEdit):
                absent = not widget.text()
            elif isinstance(widget, QLineEdit):
                absent = not widget.text().strip()
            elif isinstance(widget, QComboBox):
                absent = not widget.currentText().strip()
            else:
                absent = False
            if absent:
                missing.append(field.label)
        return missing

    @staticmethod
    def _widget(field: FieldSpec, on_change: Callable[..., None]) -> QWidget:
        kind = field.kind.lower()
        if field.choices:
            widget = QComboBox()
            widget.addItems(field.choices)
            if field.default is not None:
                widget.setCurrentText(str(field.default))
            widget.currentTextChanged.connect(on_change)
        elif kind in {"bool", "boolean", "checkbox"}:
            widget = QCheckBox()
            widget.setChecked(bool(field.default))
            widget.toggled.connect(on_change)
        elif kind in {"file", "path", "directory", "dir", "output"}:
            widget = PathEdit(kind)
            if field.default is not None:
                widget.setText(str(field.default))
            widget.changed.connect(on_change)
        else:
            widget = QLineEdit()
            if kind in {"integer", "int"}:
                widget.setValidator(QIntValidator(widget))
            elif kind in {"float", "number"}:
                widget.setValidator(QDoubleValidator(widget))
            if field.default is not None:
                widget.setText(str(field.default))
            widget.textChanged.connect(on_change)
        return widget
