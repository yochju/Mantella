#include "draggabletabwidget.hpp"

DraggableTabWidget::DraggableTabWidget(QWidget *parent) :
    QTabWidget(parent)
{
    QTabBar* tb = new DraggableTabBar;
    setTabBar(tb);
    setTabsClosable(true);
    setAcceptDrops(true);

    //setup shortcuts
    //setup remaining keys
    QShortcut *shortcut = new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_F4),this);
    connect(shortcut,SIGNAL(activated()),this,SLOT(closeTab()));
}

void DraggableTabWidget::closeTab()
{
    qDebug() << "closing tab";
    if(currentIndex() == -1) {
        return;
    }
    if(OptionsTab* optionsTab = dynamic_cast<OptionsTab*>(widget(currentIndex()))) {
        QMessageBox::information(this,"Trying to close Optionstab","You are trying to close the Optionstab. This isn't possible!");
        return;
    }
    removeTab(currentIndex());
    if(currentIndex() == -1) {
        window()->close();
    }
}

void DraggableTabWidget::dragEnterEvent(QDragEnterEvent *event)
{
    qDebug() << "widget dragenter event";
    const QMimeData *mimeData = event->mimeData();
    if(mimeData->formats().contains("action") && (mimeData->data("action") == "tab-dragging")) {
        event->accept();
    }
}

void DraggableTabWidget::dragLeaveEvent(QDragLeaveEvent *event)
{
    qDebug() << "widget dragleave event";
    event->accept();
}

void DraggableTabWidget::dragMoveEvent(QDragMoveEvent *event)
{
    qDebug() << "widget dragmove event";
    const QMimeData *mimeData = event->mimeData();
    if(mimeData->formats().contains("action") && (mimeData->data("action") == "tab-dragging")) {
        event->accept();
    }
}

void DraggableTabWidget::dropEvent(QDropEvent *event)
{
    qDebug() << "widget drop event";

    const QMimeData *mimeData = event->mimeData();
    if(mimeData->formats().contains("action") && (mimeData->data("action") == "tab-dragging")) {
        MainWindow::dndHandled = true;
        //dont need to do anything if we are dropping on the same tabbar again
        if(event->source()==tabBar()) {
            return;
        }
        //dropping on a different tabbar!
        else {
            DraggableTabBar *sourceBar = (DraggableTabBar *)event->source();
            DraggableTabWidget *sourceTabWidget = (DraggableTabWidget *)sourceBar->parentWidget();
            QWidget *draggedTab = sourceTabWidget->widget(sourceBar->currentIndex());
            QIcon icon = sourceTabWidget->tabIcon(sourceBar->currentIndex());
            QString text = sourceTabWidget->tabText(sourceBar->currentIndex());
            DraggableTabWidget *targetTabWidget = this;
            sourceTabWidget->removeTab(sourceBar->currentIndex());
            targetTabWidget->addTab(draggedTab,icon,text);
        }
        event->accept();
    }
}

void DraggableTabWidget::mouseMoveEvent(QMouseEvent *event)
{
    qDebug() << "widget mousemove event";
}

void DraggableTabWidget::mouseReleaseEvent(QMouseEvent *event)
{
    qDebug() << "widget mouserelease event";
}

void DraggableTabWidget::mousePressEvent(QMouseEvent *event)
{
    qDebug() << "widget mousepress event";
}
