#include <QtCore/QCoreApplication>

#include "DatSeg.h"
#include "In1Known.h"
#include "In1Observed.h"
#include "QAdjustFileParse.h"
int main(int argc, char* argv[])
{
    QCoreApplication a(argc, argv);
    QFile file("C:/Users/18435/Desktop/1.CPI611-CPI610往返.in1");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qDebug() << "无法打开文件：" << file.errorString();
        return 0;
    }
    QTextStream in(&file);

    QAdjustFileParse::In1::ParseIn12Enity(in);

    return a.exec();
}
