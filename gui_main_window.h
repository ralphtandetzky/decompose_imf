#pragma once

#include <QWidget>
#include <memory>

namespace gui {

namespace Ui {
class MainWindow;
}

class MainWindow : public QWidget
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private:
    std::unique_ptr<Ui::MainWindow> ui;
};


} // namespace gui
