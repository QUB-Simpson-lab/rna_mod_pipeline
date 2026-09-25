from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QImageReader, QPainter, QPixmap
from PySide6.QtWidgets import (
    QGraphicsScene, QGraphicsView, QHBoxLayout, QLabel, QPushButton,
    QSizePolicy, QToolButton, QVBoxLayout, QWidget,
)


class _ImageView(QGraphicsView):
    resized = Signal()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self.resized.emit()


class ImagePreview(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.path: Path | None = None
        self.fit_mode = True
        self.scene = QGraphicsScene(self)
        self.view = _ImageView(self.scene)
        self.view.setMinimumSize(120, 100)
        self.view.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        self.view.setBackgroundBrush(Qt.GlobalColor.white)
        self.view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.view.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.view.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.caption = QLabel()
        self.caption.setWordWrap(True)
        self.caption.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.zoom_label = QLabel()
        self.fit_button = QPushButton("Fit")
        self.actual_button = QPushButton("100%")
        self.minus_button = QToolButton()
        self.minus_button.setText("−")
        self.plus_button = QToolButton()
        self.plus_button.setText("+")
        self.open_button = QPushButton("Open image")
        self.fit_button.setToolTip("Fit the entire image without stretching it")
        self.actual_button.setToolTip("One image pixel per logical screen pixel; drag or scroll to explore")
        self.minus_button.setToolTip("Zoom out")
        self.plus_button.setToolTip("Zoom in")
        self.open_button.setToolTip("Open the original file in your system image viewer")
        self.controls = (self.fit_button, self.actual_button, self.minus_button,
                         self.plus_button, self.open_button)
        toolbar = QHBoxLayout()
        for button in self.controls[:-1]:
            toolbar.addWidget(button)
        toolbar.addWidget(self.zoom_label)
        toolbar.addStretch()
        toolbar.addWidget(self.open_button)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.caption)
        layout.addLayout(toolbar)
        layout.addWidget(self.view, 1)
        self.fit_button.clicked.connect(self.fit_image)
        self.actual_button.clicked.connect(lambda: self.set_zoom(1.0))
        self.minus_button.clicked.connect(lambda: self.set_zoom(self.scale_factor / 1.25))
        self.plus_button.clicked.connect(lambda: self.set_zoom(self.scale_factor * 1.25))
        self.open_button.clicked.connect(self._open_image)
        self.view.resized.connect(self._refit)
        self.clear()

    @property
    def scale_factor(self) -> float:
        return self.view.transform().m11()

    def clear(self) -> None:
        self.path = None
        self.fit_mode = True
        self.scene.clear()
        self.scene.setSceneRect(0, 0, 0, 0)
        self.view.resetTransform()
        self.caption.setText("Select an image to preview")
        self.caption.setToolTip("")
        self.zoom_label.clear()
        for button in self.controls:
            button.setEnabled(False)

    def set_image(self, path: Path) -> None:
        self.clear()
        reader = QImageReader(str(path))
        reader.setAutoTransform(True)
        image = reader.read()
        if image.isNull():
            self.caption.setText(f"Cannot preview {path.name}: {reader.errorString()}")
            return
        pixmap = QPixmap.fromImage(image)
        pixmap.setDevicePixelRatio(1.0)
        self.scene.addPixmap(pixmap)
        self.scene.setSceneRect(0, 0, pixmap.width(), pixmap.height())
        self.path = path
        self.caption.setText(f"{path.name}  ·  {pixmap.width()} × {pixmap.height()} px")
        self.caption.setToolTip(str(path))
        for button in self.controls:
            button.setEnabled(True)
        self.fit_image()

    def fit_image(self) -> None:
        self.fit_mode = True
        self.view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._refit()

    def _refit(self) -> None:
        if self.path is None or not self.fit_mode:
            return
        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
        self.zoom_label.setText(f"{self.scale_factor:.0%}")

    def set_zoom(self, factor: float) -> None:
        if self.path is None:
            return
        self.fit_mode = False
        self.view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.view.resetTransform()
        self.view.scale(max(0.01, min(factor, 4.0)), max(0.01, min(factor, 4.0)))
        self.zoom_label.setText(f"{self.scale_factor:.0%}")

    def _open_image(self) -> None:
        if self.path is not None:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.path)))
