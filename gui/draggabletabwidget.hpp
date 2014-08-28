#ifndef DRAGGABLETABWIDGET_HPP
#define DRAGGABLETABWIDGET_HPP

#include <QDebug>

#include <QTabWidget>
#include <QtGui>
#include <QShortcut>
#include <QMessageBox>
#include "draggabletabbar.hpp"

class DraggableTabWidget : public QTabWidget
{
    Q_OBJECT
public:
    explicit DraggableTabWidget(QWidget *parent = 0);

signals:

public slots:
    //shortcuts
    void closeTab();

protected:
    virtual void dragEnterEvent(QDragEnterEvent *event);
    virtual void dragLeaveEvent(QDragLeaveEvent *event);
    virtual void dragMoveEvent(QDragMoveEvent *event);
    virtual void dropEvent(QDropEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
};

#endif // DRAGGABLETABWIDGET_HPP
