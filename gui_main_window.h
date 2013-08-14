#pragma once

#include <QWidget>
#include <memory>

namespace gui {

class MainWindow : public QWidget
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void optimize();
    
private:
    struct Impl;
    std::unique_ptr<Impl> m;
};


} // namespace gui
