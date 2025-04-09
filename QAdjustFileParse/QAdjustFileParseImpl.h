#pragma once
#include "QAdjustFileParse.h"
namespace QAdjustFileParse {
/**
 * @brief 空格替换分割
 */
QStringList SpaceReplacSplit(QString input);

/**
 * @brief 两数绝对值的均值
 * @param num1
 * @param num2
 * @return 结果
 */
double absoluteMean(double num1, double num2);
namespace GSI {

    /**
     * @brief GSI数据解析
     * @param stream 文件流
     * @return QList{<QList{<QList{{目标，方向，高差，距离}},<目标，高差之差，高差之差累计差，前后视距累计差，距起始点距离，高程>>},<起始点，起始点高程>>}
     */
    std::optional<QList<QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>>>> ParseGSI(QTextStream& stream);

    /**
     * @brief GSI文件测回数据解析
     * @param stream 文件流
     * @return <QList{{目标，方向，高差，距离}},<目标，高差之差，高差之差累计差，前后视距累计差，距起始点距离，高程>>
     */
    std::optional<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>> ParseGSIRound(QTextStream& stream);

    /**
     * @brief GSI文件点名字段处理
     * @param input
     * @return 点名
     */
    QString ParseGSIPName(QString input);

    /**
     * @brief GSI文件数据字段处理
     * @param input
     * @return 数值
     */
    double ParseGSIValue(QString input);
}

namespace Dat {

    /**
     * @brief dat数据解析
     * @param stream 文件流
     * <QList{{目标，方向，高差，距离}},<目标，高程>>--->测回数据
     * @return QList<<QList{{目标，方向，高差，距离}},<目标，高程>>,<目标，高程>>--->测段数据
     */
    std::optional<QList<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>>> ParseDat(QTextStream& stream);

    /**
     * @brief dat测段数据解析
     * @param stream 文件流
     * @return  <QList{<QList{{目标，方向，高差，距离}},<目标，高程>>},<起始点，起始点高程>>
     */
    std::optional<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>> ParseDatSeg(QTextStream& stream);

    /**
     * @brief dat文件测回数据解析
     * @param stream 文件流
     * @return <QList{{目标，方向，高差，距离}},<目标，高程>>
     */
    std::optional<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>> ParseDatRound(QTextStream& stream);

    /**
     * @brief dat文件有效行解析
     * @param line 行数据
     * @return {目标，方向，高差，距离}
     */
    std::optional<std::tuple<QString, QString, double, double>> ParseDatRoundLine(QString line);

    /**
     * @brief dat文件高程行解析
     * @param line 行数据
     * @return <目标，高程>
     */
    std::optional<QPair<QString, double>> ParseDatRoundAltitudeLine(QString line);

    /**
     * @brief 是否是dat文件的有效行
     * @param line 行数据
     * @return
     */
    bool IsValidDatLine(QString line);
}

namespace SUC {
    /**
     * @brief suc文件测回数据解析
     * @param stream 文件流
     * @return QList<QPair<点名，{QPair<{水平角(左)，天顶距1(左)，{斜距1(左)},{水平角(右)，天顶距1(右)，斜距1(右)}>,仪器高，棱镜高}>>
     */
    std::optional<QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>> ParseSucRound(QTextStream& stream);

    std::optional<QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>> ParseSucRoundEveryOrient(QTextStream& stream);

    /**
     * @brief suc文件测回数据解析
     * @param stream 文件流
     * @return QPair<测站名,QList<QPair<点名，{QPair<{水平角(左)，天顶距1(左)，{斜距1(左)},{水平角(右)，天顶距1(右)，斜距1(右)}>,仪器高，棱镜高}>>>
     */
    std::optional<QPair<QString, QList<QPair<int, QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>>>>> ParseSuc(QTextStream& stream);
}
namespace In2 {
    std::optional<QPair<QString, QList<std::tuple<QString, QString, double>>>> ParseIn2StnData(QTextStream& stream);
}
namespace In1 {

    /**
     * @brief int文件解析
     * @param stream 文件流
     * @return QMap{<起始点，目标点>，<高差，距离>}
     */
    std::optional<QMap<QPair<QString, QString>, QPair<double, double>>> ParseIn1(QTextStream& stream);
}
};
