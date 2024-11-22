#pragma once

#include "qadjustfileparse_global.h"

#include <QFile>
#include <QMap>
#include <QObject>
#include <QPair>
#include <QTextStream>

#include <Eigen/Dense>

#include "DatSeg.h"
#include "GSISeg.h"
#include "In1Observed.h"
#include "In2Observed.h"
#include "SucPoint.h"

namespace QAdjustFileParse {

namespace Dat {

    /**
     * @brief dat数据解析
     * @param stream 文件流
     * <QList{{目标，方向，高差，距离}},<目标，高程>>--->测回数据
     * @return QList<<QList{{目标，方向，高差，距离}},<目标，高程>>,<目标，高程>>--->测段数据
     */
    QADJUSTFILEPARSE_EXPORT QList<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>> ParseDat(QTextStream& stream);

    /**
     * @brief 解析dat文件为对象
     * @param stream 文件流
     * @return
     */
    QADJUSTFILEPARSE_EXPORT QList<DatSeg> ParseDat2Entity(QTextStream& stream);
}

namespace GSI {
    /**
     * @brief GSI数据解析
     * @param stream 文件流
     * @return QList{<QList{<QList{{目标，方向，高差，距离}},<目标，高差之差，高差之差累计差，前后视距累计差，距起始点距离，高程>>},<起始点，起始点高程>>}
     */
    QADJUSTFILEPARSE_EXPORT QList<QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>>> ParseGSI(QTextStream& stream);

    /**
     * @brief GSI数据解析为对象
     * @param stream
     * @return
     */
    QADJUSTFILEPARSE_EXPORT QList<GSISeg> ParseGSI2Entity(QTextStream& stream);
}
namespace In1 {
    /**
     * @brief int文件解析
     * @param stream 文件流
     * @return QMap{<起始点，目标点>，<高差，距离>}
     */
    QADJUSTFILEPARSE_EXPORT QMap<QPair<QString, QString>, QPair<double, double>> ParseIn1(QTextStream& stream);

    /**
     * @brief in1解析为往返测数据
     * @param stream 文件流
     * @return  QMap{<起始点，目标点>，<<往测高差，往测距离>,<返测高差，返测距离>>}
     */
    QADJUSTFILEPARSE_EXPORT QMap<QPair<QString, QString>, QPair<QPair<double, double>, QPair<double, double>>> ParseIn1EveryOrient(QTextStream& stream);

    QADJUSTFILEPARSE_EXPORT QList<In1Observed> ParseIn12Enity(QTextStream& stream);
}

namespace SUC {
    /**
     * @brief suc文件测回数据解析
     * @param stream 文件流
     * @return QPair<测站名,QList<QPair<点名，{QPair<{水平角(左)，天顶距1(左)，{斜距1(左)},{水平角(右)，天顶距1(右)，斜距1(右)}>,仪器高，棱镜高}>>>
     */
    QADJUSTFILEPARSE_EXPORT QPair<QString, QList<QPair<int, QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>>>> ParseSuc(QTextStream& stream);

    /**
     * @brief suc文件测回数据解析
     * @param stream 文件流
     * @return QMap<测站名,QPair<测回，测回数据>>
     */
    QADJUSTFILEPARSE_EXPORT QMap<QString, QList<QPair<int, QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>>>> ParseSucEveryOrient(QTextStream& stream);

    /**
     * @brief suc数据解析为对象
     * @param stream
     * @return
     */
    QADJUSTFILEPARSE_EXPORT QList<SucPoint> ParseSuc2Entity(QTextStream& stream);
}

namespace In2 {
    /**
     * @brief in2文件解析
     * @param stream 文件流
     * @return tuple(tuple(方向中误差，方位角误差，测边固定误差，比例误差)}，QMap{测站号，Eigen::Vector2d<x,y>},QList<QPair<测站号，QList<tuple(观测名，类型，观测值)>>>)
     */
    QADJUSTFILEPARSE_EXPORT QList<QPair<QString, QList<std::tuple<QString, QString, double>>>> ParseIn2(QTextStream& stream);

    /**
     * @brief in2解析为对象
     * @param stream
     * @return
     */
    QADJUSTFILEPARSE_EXPORT QList<QPair<QString, QList<In2Observed>>> ParseIn22Entity(QTextStream& stream);
}
namespace Gra {

    QADJUSTFILEPARSE_EXPORT QPair<QMap<QString, QPair<Eigen::Vector2d, bool>>, QList<QPair<QString, QString>>> ParseGra(QTextStream& stream);
}
};
