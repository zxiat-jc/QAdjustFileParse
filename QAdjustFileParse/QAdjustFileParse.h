#pragma once

#include "qadjustfileparse_global.h"

#include <QFile>
#include <QMap>
#include <QObject>
#include <QPair>
#include <QTextStream>

#include <Eigen/Dense>

#include "DatSeg.h"
#include "GSIRound.h"
#include "GSISeg.h"
#include "In1Observed.h"
#include "In2Observed.h"
#include "SucPoint.h"

namespace QAdjustFileParse {

namespace Dat {

    /**
     * @brief 解析dat文件为对象
     * @param stream 文件流
     * @return
     */
    QADJUSTFILEPARSE_EXPORT std::optional<QList<DatSeg>> ParseDat2Entity(QTextStream& stream);
}

namespace GSI {

    /**
     * @brief 解析gsi文件为对象
     * @param stream 文件流
     * @return
     */
    QADJUSTFILEPARSE_EXPORT std::optional<QList<GSISeg>> ParseGSI2Entity(QTextStream& stream);

}
namespace In1 {

    /**
     * @brief in1解析为往返测数据
     * @param stream 文件流
     * @return  QMap{<起始点，目标点>，<<往测高差，往测距离>,<返测高差，返测距离>>}
     */
    QADJUSTFILEPARSE_EXPORT std::optional<QMap<QPair<QString, QString>, QPair<QPair<double, double>, QPair<double, double>>>> ParseIn1EveryOrient(QTextStream& stream);

    QADJUSTFILEPARSE_EXPORT std::optional<QList<In1Observed>> ParseIn12Entity(QTextStream& stream);
}

namespace SUC {
    QADJUSTFILEPARSE_EXPORT std::optional<QList<SucPoint>> ParseSuc2Entity(QTextStream& stream);

    /**
     * @brief suc文件测回数据解析
     * @param stream 文件流
     * @return QMap<测站名,QPair<测回，测回数据>>
     */
    QADJUSTFILEPARSE_EXPORT std::optional<QMap<QString, QList<QPair<int, QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>>>>> ParseSucEveryOrient(QTextStream& stream);
}

namespace In2 {
    /**
     * @brief in2文件解析
     * @param stream 文件流
     * @return tuple(tuple(方向中误差，方位角误差，测边固定误差，比例误差)}，QMap{测站号，Eigen::Vector2d<x,y>},QList<QPair<测站号，QList<tuple(观测名，类型，观测值)>>>)
     */
    QADJUSTFILEPARSE_EXPORT std::optional<QList<QPair<QString, QList<std::tuple<QString, QString, double>>>>> ParseIn2(QTextStream& stream);

    /**
     * @brief in2解析为对象
     * @param stream
     * @return
     */
    QADJUSTFILEPARSE_EXPORT std::optional<QList<QPair<QString, QList<In2Observed>>>> ParseIn22Entity(QTextStream& stream);
}
namespace Gra {
    QADJUSTFILEPARSE_EXPORT QPair<QMap<QString, QPair<Eigen::Vector2d, bool>>, QList<QPair<QString, QString>>> ParseGra(QTextStream& stream);
}
}
