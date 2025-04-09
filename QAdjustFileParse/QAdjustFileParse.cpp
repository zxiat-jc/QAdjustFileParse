#include "QAdjustFileParse.h"
#include "QAdjustFileParseImpl.h"

#include <QDebug>
#include <QFile>
#include <QJsonArray>
#include <QMap>
#include <QPair>
#include <QTextStream>

#include <regex>

#include "ConfigWrap.h"
#include "DatRound.h"
#include "QtUtils.h"

using namespace QAdjustFileParse;

std::optional<QList<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>>> QAdjustFileParse::Dat::ParseDat(QTextStream& in)
{
    QList<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>> dat;
    while (!in.atEnd()) {
        auto&& opt = Dat::ParseDatSeg(in);
        if (!opt.has_value()) {
            return std::nullopt;
        }
        dat.push_back(opt.value());
    }
    return dat;
}

 std::optional<QList<DatSeg>> QAdjustFileParse::Dat::ParseDat2Entity(QTextStream& stream)
{
    QList<DatSeg> segs;
    auto&& opt = QAdjustFileParse::Dat::ParseDat(stream);
    if (!opt.has_value()) {
        return std::nullopt;
    }
    auto&& dat = opt.value();
    for (int i = 0; i < dat.size(); i++) {
        // 测段
        auto&& datSeg = dat[i];
        // 起始点
        QString startPoint = datSeg.second.first;
        // 起始点高程
        double startAltitude = datSeg.second.second;
        QList<DatRound> rounds;
        // 结束点
        QString endPoint;
        // 高差，距离
        double heightDiff = 0, distance = 0;
        // 处理测段内的测回
        for (int j = 0; j < datSeg.first.size(); j++) {
            auto&& datRound = datSeg.first[j];
            DatRound round(datRound);
            endPoint = round.fName();
            distance += round.distance();
            heightDiff += round.heightDiff();
            rounds.append(round);
        }
        DatSeg seg(rounds, startAltitude, startPoint, endPoint, heightDiff, distance);
        segs.append(seg);
    }
    return segs;
}

std::optional<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>> QAdjustFileParse::Dat::ParseDatSeg(QTextStream& in)
{
    // 测段 数据结构
    QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>> seg;
    // 测段起始点数据
    QPair<QString, double> start;
    // 判断时是否为新测段
    bool s = false;
    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine();
        if (!line.startsWith("For M5")) {
            qDebug() << "DAT测回数据有误" << line;
            return std::nullopt;
        }
        if (line.contains("Start-Line")) {
            s = true;
            auto&& opt = Dat::ParseDatRoundAltitudeLine(in.readLine());
            if (!opt.has_value()) {
                return std::nullopt;
            }
            start = opt.value();
        } else if (line.contains("End-Line")) {
            s = false;
            break;
        } else {
            if (!Dat::IsValidDatLine(line)) {
                continue;
            }
            in.seek(pos);
        }

        if (s) {
            auto&& opt = Dat::ParseDatRound(in);
            if (!opt.has_value()) {
                return std::nullopt;
            }
            auto&& round = opt.value();
            seg.push_back(round);
        }
    }
    return qMakePair(seg, start);
}

std::optional<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>> QAdjustFileParse::Dat::ParseDatRound(QTextStream& in)
{
    QList<std::tuple<QString, QString, double, double>> roundData;
    QPair<QString, double> targetHeight;
    // 是否是测回起始点所在行
    std::optional<QString> startPoint { std::nullopt };
    while (!in.atEnd()) {
        QString line = in.readLine();
        if (!Dat::IsValidDatLine(line)) {
            continue;
        }
        auto&& opt = Dat::ParseDatRoundLine(line);
        if (!opt.has_value()) {
            return std::nullopt;
        }
        auto&& lineData = opt.value();

        if (!startPoint.has_value()) {
            startPoint = std::get<0>(lineData);
            roundData.push_back(lineData);
        } else if (startPoint == std::get<0>(lineData)) {
            roundData.push_back(lineData);
            auto&& opt = ParseDatRoundAltitudeLine(in.readLine());
            if (!opt.has_value()) {
                return std::nullopt;
            }
            targetHeight = opt.value();
            break;
        } else {
            roundData.push_back(lineData);
        }
    }
    return qMakePair(roundData, targetHeight);
}

std::optional<std::tuple<QString, QString, double, double>> QAdjustFileParse::Dat::ParseDatRoundLine(QString line)
{
    QStringList parts = line.split("|");
    if (parts.size() != 7) {
        qDebug() << "DAT测回数据有误" << line;
        return std::nullopt;
    }
    parts.removeFirst();
    parts.removeFirst();
    QString target, orient;
    double heightDiff, distance;

    auto&& targetSection = SpaceReplacSplit(parts[0]);
    auto&& heightSection = SpaceReplacSplit(parts[1]);
    auto&& distanceSection = SpaceReplacSplit(parts[2]);
    if (!(targetSection.size() >= 2) || heightSection.size() != 3 || distanceSection.size() != 3) {
        qDebug() << "DAT测回数据有误" << line;
        return std::nullopt;
    }
    target = targetSection[1];
    orient = heightSection[0];
    heightDiff = heightSection[1].toDouble();
    distance = distanceSection[1].toDouble();
    return std::make_tuple(target, orient, heightDiff, distance);
}

std::optional<QPair<QString, double>> QAdjustFileParse::Dat::ParseDatRoundAltitudeLine(QString line)
{
    QStringList parts = line.split("|");
    if (parts.size() != 7) {
        qDebug() << "DAT测回高程行数据有误" << line;
        return std::nullopt;
    }
    parts.removeFirst();
    parts.removeFirst();

    auto&& targetSection = SpaceReplacSplit(parts[0]);
    auto&& altitudeSection = SpaceReplacSplit(parts[3]);
    if (targetSection.size() < 2 || altitudeSection.size() < 2) {
        qDebug() << "DAT测回高程行数据有误" << line;
        return std::nullopt;
    }
    QString target = targetSection[1];
    double altitude = altitudeSection[1].toDouble();
    return qMakePair(target, altitude);
}

std::optional<QMap<QPair<QString, QString>, QPair<double, double>>> QAdjustFileParse::In1::ParseIn1(QTextStream& in)
{
    // 是否读到未知点所在行
    bool s = false;
    // 已知点数据
    // QMap<QString, double> knownPoints;
    // 来回测数据、单向数据
    QMap<QPair<QString, QString>, QPair<double, double>> data;
    int lineNum = 0;
    while (!in.atEnd()) {
        lineNum++;
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) {
            continue;
        }
        QStringList splitString = line.split(",");
        if (splitString.size() == 4) {
            s = true;
        }

        if (s) {
            if (splitString.size() == 4) {
                QString startp = splitString[0];
                QString endp = splitString[1];
                double heightDiff = splitString[2].toDouble();
                double distance = splitString[3].toDouble();

                // 检查之前存储数据 是否与当前数据配对
                if (data.contains(qMakePair(endp, startp))) {
                    // 处理后数据
                    double fHeightDiff = data[qMakePair(endp, startp)].first;
                    double fDistance = data[qMakePair(endp, startp)].second;
                    data[qMakePair(endp, startp)].first = fHeightDiff < 0 ? -(absoluteMean(fHeightDiff, heightDiff)) : absoluteMean(fHeightDiff, heightDiff);
                    data[qMakePair(endp, startp)].second = absoluteMean(fDistance, distance);
                } else {
                    data[qMakePair(startp, endp)].first = heightDiff;
                    data[qMakePair(startp, endp)].second = distance;
                }
            } else {
                qDebug() << "In1数据有误" << line;
                return std::nullopt;
            }
        }
    }
    return data;
}

std::optional<QMap<QPair<QString, QString>, QPair<QPair<double, double>, QPair<double, double>>>> QAdjustFileParse::In1::ParseIn1EveryOrient(QTextStream& in)
{
    // 是否读到未知点所在行
    bool s = false;
    // 来回测数据、单向数据
    QMap<QPair<QString, QString>, QPair<QPair<double, double>, QPair<double, double>>> data;
    int lineNum = 0;
    while (!in.atEnd()) {
        lineNum++;
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) {
            continue;
        }
        QStringList splitString = line.split(",");

        if (splitString.size() == 4) {
            s = true;
        }
        if (s) {
            if (splitString.size() != 4) {
                return std::nullopt;
            }
            QString startp = splitString[0];
            QString endp = splitString[1];
            double heightDiff = splitString[2].toDouble();
            double distance = splitString[3].toDouble();

            // 检查之前存储数据 是否与当前数据配对
            if (data.contains(qMakePair(endp, startp))) {
                data[qMakePair(endp, startp)].second = qMakePair(heightDiff, distance);
            } else {
                data[qMakePair(startp, endp)].first = qMakePair(heightDiff, distance);
            }
        }
    }
    // 移除未匹配数据
    for (auto&& pair : data.keys()) {
        if (data[pair].second == QPair<double, double>()) {
            data.remove(pair);
        }
    }

    return data;
}

 std::optional<QList<In1Observed>> QAdjustFileParse::In1::ParseIn12Entity(QTextStream& stream)
{
    QList<In1Observed> obss;
    auto&& opt = QAdjustFileParse::In1::ParseIn1(stream);
    if (!opt.has_value()) {
        return std::nullopt;
    }
    auto&& in1 = opt.value();
    auto&& observedDatas = in1;
    for (auto&& pointPair : observedDatas.keys()) {
        auto&& start = pointPair.first;
        auto&& end = pointPair.second;
        auto&& hDiff = observedDatas[pointPair].first;
        auto&& distance = observedDatas[pointPair].second;
        In1Observed in1Ob(start, end, hDiff, distance);
        obss.append(in1Ob);
    }
    return obss;
}

bool QAdjustFileParse::Dat::IsValidDatLine(QString line)
{
    if (line.contains("####") || line.contains("Measurement repeated") || line.contains("Station repeated") || line.contains("Reading E327") || line.contains("Sh") || line.contains("Db")) {
        return false;
    }
    return true;
}

QStringList QAdjustFileParse::SpaceReplacSplit(QString input)
{
    std::regex pattern("\\s+");
    std::string result = std::regex_replace(input.toStdString(), pattern, " ");
    return QString::fromStdString(result).trimmed().split(" ");
}

double QAdjustFileParse::absoluteMean(double num1, double num2)
{
    return (abs(num1) + abs(num2)) / 2;
}

std::optional<QList<QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>>>> QAdjustFileParse::GSI::ParseGSI(QTextStream& in)
{
    // 存储不同测段的已知点高程
    QMap<QString, double> knownHeight;
    // 存储多个测段信息
    QList<QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>>> gsis;
    while (!in.atEnd()) {
        // 测段已知点高程数据
        QPair<QString, double> start;
        // 存储多个测回信息
        QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>> rounds;
        // 是否开始测段测量
        bool s = false;
        while (!in.atEnd()) {
            auto&& pos = in.pos();
            QString line = in.readLine();
            if (!line.startsWith("*")) {
                qDebug() << "GSI数据有误" << line;
                return std::nullopt;
            }

            QStringList list = SpaceReplacSplit(line);
            if (s) {
                in.seek(pos);
                if (list.size() == 1) {
                    break;
                } else {
                    auto&& opt = GSI::ParseGSIRound(in);
                    if (!opt.has_value()) {
                        return std::nullopt;
                    }
                    auto&& round = opt.value();
                    rounds.append(round);
                    continue;
                }
            }
            // 读到测段起始点
            if (list.size() == 1) {
                s = true;
                // 处理起始点数据
                QString startLine = in.readLine();
                QStringList startList = SpaceReplacSplit(startLine);
                if (startList.size() == 2) {
                    start.first = GSI::ParseGSIPName(startList[0]);
                    start.second = GSI::ParseGSIValue(startList[1]);
                    knownHeight[start.first] = start.second;
                } else if (startList.size() == 6) {
                    auto&& pname = GSI::ParseGSIPName(startList[0]);
                    if (knownHeight.contains(pname)) {
                        start.first = pname;
                        start.second = knownHeight[startList[0]];
                    }
                } else {
                    qDebug() << "GSI测段起始点数据有误" << startLine;
                    return std::nullopt;
                }
            }
        }
        gsis.push_back(qMakePair(rounds, start));
    }
    return gsis;
}

std::optional<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>> QAdjustFileParse::GSI::ParseGSIRound(QTextStream& in)
{
    QList<std::tuple<QString, bool, double, double>> rdatas;
    std::tuple<QString, double, double, double, double, double> sdata;

    QStringList sLine;
    QList<QStringList> rLines;
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList lineList = SpaceReplacSplit(line);
        if (lineList.size() == 5) {
            rLines.append(lineList);
        } else if (lineList.size() == 6) {
            sLine = lineList;
            break;
        } else {
            return std::nullopt;
        }
    }
    // 目标
    QString targetPName = GSI::ParseGSIPName(sLine[0]);
    // 高差之差
    double heightDiffDiff = GSI::ParseGSIValue(sLine[1]);
    // 高差之差累计差
    double heightDiffDiffAcc = GSI::ParseGSIValue(sLine[2]);
    // 前后视距累计差
    double distAccDiff = GSI::ParseGSIValue(sLine[3]);
    // 距起始点距离
    double distFromStart = GSI::ParseGSIValue(sLine[4]);
    // 高程
    double altitude = GSI::ParseGSIValue(sLine[5]);
    sdata = std::make_tuple(targetPName, heightDiffDiff, heightDiffDiffAcc, distAccDiff, distFromStart, altitude);

    // 是否有前侧点数据
    bool s = false;
    // 测回四条数据 是否与目标点匹配
    for (auto&& rLine : rLines) {
        // 前进目标方向为true
        bool orient = false;
        QString target = GSI::ParseGSIPName(rLine[0]);
        // 标尺读数
        double rulerNum = GSI::ParseGSIValue(rLine[1]);
        double dist = GSI::ParseGSIValue(rLine[2]);
        targetPName == target ? orient = true, s = true : orient = false;
        rdatas.append(std::make_tuple(target, orient, rulerNum, dist));
    }

    // 未读到四个数据
    if (rdatas.size() != 4 || !s) {
        qDebug() << "GSI测回数据有误";
        qDebug() << rLines;
        return std::nullopt;
    }
    // 四条数据不匹配
    if (std::get<0>(rdatas[0]) != std::get<0>(rdatas[3]) || std::get<0>(rdatas[1]) != std::get<0>(rdatas[2])) {
        qDebug() << "GSI测回数据有误";
        qDebug() << rLines;
        return std::nullopt;
    }

    return qMakePair(rdatas, sdata);
}

QString QAdjustFileParse::GSI::ParseGSIPName(QString input)
{
    QString right = input.right(8);
    /* while (right.startsWith("0")) {
         right = right.right(right.size() - 1);
     }*/
    return right;
}

double QAdjustFileParse::GSI::ParseGSIValue(QString input)
{
    double data = input.right(8).toInt() / 100000.0;
    if (input.contains("-")) {
        return -data;
    }
    return data;
}

std::optional<QPair<QString, QList<QPair<int, QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>>>>> QAdjustFileParse::SUC::ParseSuc(QTextStream& in)
{
    // 测站名
    QString deviceName;
    // 是否读到测回行
    bool s = false;
    QList<QPair<int, QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>>> rounds;
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList lineList = line.split(",");
        if (!s) {
            if (lineList.size() == 4) {
                deviceName = lineList[0].trimmed();
            }
            if (lineList.size() == 1) {
                s = true;
            }
            continue;
        }
        if (s) {
            int surveyTime = lineList[0].trimmed().toInt();
            auto&& opt = SUC::ParseSucRound(in);
            if (!opt.has_value()) {
                qDebug() << "SUC数据有误"
                         << "测站名：" << deviceName << "测回：" << surveyTime;
                return std::nullopt;
            }
            auto&& round = opt.value();
            rounds.append(qMakePair(surveyTime, round));
        }
    }
    return qMakePair(deviceName, rounds);
}

 std::optional<QList<SucPoint>> QAdjustFileParse::SUC::ParseSuc2Entity(QTextStream& stream)
{
    QList<SucPoint> points;
    auto&& opt = QAdjustFileParse::SUC::ParseSuc(stream);
    if (!opt.has_value()) {
        return std::nullopt;
    }
    auto&& suc = opt.value();
    QString stnName = suc.first;
    for (auto&& round : suc.second) {
        int surveyTime = round.first;
        for (auto&& point : round.second) {
            QString pname = point.first;
            auto&& [hf, vf, sf] = std::get<0>(point.second).first;
            auto&& [hr, vr, sr] = std::get<0>(point.second).second;
            auto&& deviceHeight = std::get<1>(point.second);
            auto&& prismHeight = std::get<2>(point.second);
            SucPoint point(stnName, surveyTime, pname, Utils::Geometry::Fun::HV::AMS2A(hf), Utils::Geometry::Fun::HV::AMS2A(vf), sf, Utils::Geometry::Fun::HV::AMS2A(hr), Utils::Geometry::Fun::HV::AMS2A(vr), sr, deviceHeight, prismHeight);
            points.append(point);
        }
    }
    return points;
}

std::optional<QMap<QString, QList<QPair<int, QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>>>>> QAdjustFileParse::SUC::ParseSucEveryOrient(QTextStream& in)
{
    // 测站名
    QString deviceName;
    // 是否读到测回行
    bool s = false;
    // 测回
    int surveyTime = 0;
    // 是否是新的测站
    bool isNewStn = false;
    // 所有测回数据
    QList<QPair<int, QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>>> rounds;
    // 结果
    QMap<QString, QList<QPair<int, QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>>>> result;
    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine();
        QStringList lineList = line.split(",");
        if (!s) {
            if (lineList.size() == 4) {
                if (isNewStn) {
                    deviceName = lineList[0].trimmed();
                    isNewStn = true;
                    rounds.clear();
                } else {
                    deviceName = lineList[0].trimmed();
                    isNewStn = true;
                }
            }
            if (lineList.size() == 1) {
                s = true;
                surveyTime = lineList[0].trimmed().toInt();
            }
        }
        if (s && isNewStn) {
            auto&& opt = QAdjustFileParse::SUC::ParseSucRoundEveryOrient(in);
            if (!opt.has_value()) {
                qDebug() << "SUC数据处理有误"
                         << "测站名：" << deviceName << "测回：" << surveyTime;
                return std::nullopt;
            }
            auto&& round = opt.value();
            rounds.append(qMakePair(surveyTime, round));
            s = false;
        }
        if (in.atEnd()) {
            result[deviceName] = rounds;
            return result;
        }
    }
    return result;
}

std::optional<QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>> QAdjustFileParse::SUC::ParseSucRound(QTextStream& in)
{
    // 上半测回点集map key:点名 value:点序号 序号是leftDatas的索引
    QMap<QString, int> leftPnames;
    QMap<QString, int> rightPnames;
    int leftIndexs = 0;
    int rightIndexs = 0;
    // 上半测回数据
    QList<QPair<QString, std::tuple<std::tuple<double, double, double>, double, double>>> leftDatas;
    // 下半测回数据
    QList<QPair<QString, std::tuple<std::tuple<double, double, double>, double, double>>> rightDatas;
    // 上半测回是否结束
    bool isEnd = false;

    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine();
        if (line.contains("End")) {
            break;
        }
        QStringList lineList = line.split(",");
        if (lineList.size() == 6) {
            QString pname = lineList[0].trimmed();
            double h = lineList[1].trimmed().toDouble();
            double v = lineList[2].trimmed().toDouble();
            double s = lineList[3].trimmed().toDouble();
            double deviceHeight = lineList[4].trimmed().toDouble();
            double prismHeight = lineList[5].trimmed().toDouble();
            QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>> round;
            // 上半测回
            if (!isEnd) {
                // 当测到已存在的点时，该测回结束
                if (leftPnames.contains(pname)) {
                    auto&& [ph, pv, ps] = std::get<0>(leftDatas[leftPnames[pname]].second);
                    std::get<0>(leftDatas[leftPnames[pname]].second) = std::make_tuple((h + ph) / 2, (v + pv) / 2, (s + ps) / 2);
                    isEnd = true;
                } else {
                    leftDatas.append(qMakePair(pname, std::make_tuple(std::make_tuple(h, v, s), deviceHeight, prismHeight)));
                    leftPnames[pname] = leftIndexs;
                    leftIndexs++;
                }
            } else {
                if (rightPnames.contains(pname)) {
                    auto&& [ph, pv, ps] = std::get<0>(rightDatas[rightPnames[pname]].second);
                    std::get<0>(rightDatas[rightPnames[pname]].second) = std::make_tuple((h + ph) / 2, (v + pv) / 2, (s + ps) / 2);
                } else {
                    rightDatas.append(qMakePair(pname, std::make_tuple(std::make_tuple(h, v, s), deviceHeight, prismHeight)));
                    rightIndexs++;
                    rightPnames[pname] = rightIndexs;
                }
            }
        } else if (lineList.size() == 1) {
            in.seek(pos);
            break;
        } else {
            return std::nullopt;
        }
    }

    if (leftDatas.size() != rightDatas.size()) {
        return std::nullopt;
    }

    QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>> rounds;
    auto&& leftIter = leftDatas.begin();
    auto&& rightIter = rightDatas.begin();
    while (leftIter != leftDatas.end() && rightIter != rightDatas.end()) {
        auto&& ldata = *leftIter;
        auto&& rdata = *rightIter;

        auto&& [lh, lv, ls] = std::get<0>(ldata.second);
        auto&& [rh, rv, rs] = std::get<0>(rdata.second);

        auto&& deviceHeight = std::get<1>(ldata.second);
        auto&& prismHeight = std::get<2>(ldata.second);

        rounds.append(qMakePair(ldata.first, std::make_tuple(qMakePair(std::make_tuple(lh, lv, ls), std::make_tuple(rh, rv, rs)), deviceHeight, prismHeight)));

        leftIter++;
        rightIter++;
    }

    return rounds;
}

std::optional<QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>>> QAdjustFileParse::SUC::ParseSucRoundEveryOrient(QTextStream& in)
{
    // 上半测回点集map key:点名 value:点序号 序号是leftDatas的索引
    QMap<QString, int> leftPnames;
    QMap<QString, int> rightPnames;
    int leftIndexs = 0;
    int rightIndexs = 0;
    // 上半测回数据
    QList<QPair<QString, std::tuple<double, double, double>>> leftDatas;
    // 下半测回数据
    QList<QPair<QString, std::tuple<double, double, double>>> rightDatas;
    // 上半测回是否结束
    bool isEnd = false;

    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine();
        QStringList lineList = line.split(",");
        if (lineList.size() == 6) {
            QString pname = lineList[0].trimmed();
            double h = lineList[1].trimmed().toDouble();
            double v = lineList[2].trimmed().toDouble();
            double s = lineList[3].trimmed().toDouble();

            // 上半测回
            if (!isEnd) {
                // 当测到已存在的点时，该测回结束
                leftDatas.append(qMakePair(pname, std::make_tuple(h, v, s)));
                if (leftPnames.contains(pname)) {
                    isEnd = true;
                } else {
                    leftPnames[pname] = leftIndexs;
                    leftIndexs++;
                }
            } else {
                rightDatas.append(qMakePair(pname, std::make_tuple(h, v, s)));
                rightIndexs++;
                rightPnames[pname] = rightIndexs;
            }
        } else if (lineList.size() == 1) {
            in.seek(pos);
            break;
        }
    }

    QList<QPair<QString, QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>>> rounds;
    if (leftDatas.size() != rightDatas.size()) {
        return std::nullopt;
    }
    auto&& leftIter = leftDatas.begin();
    auto&& rightIter = rightDatas.rbegin();

    while (leftIter != leftDatas.end() && rightIter != rightDatas.rend()) {
        auto&& ldata = *leftIter;
        auto&& rdata = *rightIter;

        auto&& [lh, lv, ls] = ldata.second;
        auto&& [rh, rv, rs] = rdata.second;
        if (ldata.first == rdata.first) {
            rounds.append(qMakePair(ldata.first, qMakePair(std::make_tuple(lh, lv, ls), std::make_tuple(rh, rv, rs))));
        }
        leftIter++;
        rightIter++;
    }

    return rounds;
}

std::optional<QList<QPair<QString, QList<std::tuple<QString, QString, double>>>>> QAdjustFileParse::In2::ParseIn2(QTextStream& in)
{
    // 观测数据
    QList<QPair<QString, QList<std::tuple<QString, QString, double>>>> stnDatas;
    // 是否读到测站行
    bool s = false;
    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine().trimmed();
        // 空行跳过
        if (line == "") {
            continue;
        }
        QStringList lineList = line.split(",");
        if (!s) {
        }

        if (lineList.size() == 1) {
            in.seek(pos);
            auto&& opt = ParseIn2StnData(in);
            if (!opt.has_value()) {
                return std::nullopt;
            }
            auto&& stnData = opt.value();
            stnDatas.append(stnData);
            s = true;
        }
    }
    return stnDatas;
}

std::optional<QList<QPair<QString, QList<In2Observed>>>> QAdjustFileParse::In2::ParseIn22Entity(QTextStream& stream)
{
    auto&& opt = QAdjustFileParse::In2::ParseIn2(stream);
    if (!opt.has_value()) {
        return std::nullopt;
    }
    auto&& ObservedGroups = opt.value();
    QList<QPair<QString, QList<In2Observed>>> result;
    for (auto&& ObservedGroup : ObservedGroups) {
        QPair<QString, QList<In2Observed>> stnData;
        QList<In2Observed> obss;
        QString stn = ObservedGroup.first;
        for (auto&& obsData : ObservedGroup.second) {
            auto&& [pname, type, value] = obsData;
            auto&& obsType = S2E(In2Observed::ObservedValueType, type);
            In2Observed obs(stn, pname, S2E(In2Observed::ObservedValueType, type), value);
            obss.append(obs);
        }
        stnData = qMakePair(stn, obss);
        result.append(stnData);
    }
    return result;
}

std::optional<QPair<QString, QList<std::tuple<QString, QString, double>>>> QAdjustFileParse::In2::ParseIn2StnData(QTextStream& in)
{
    QPair<QString, QList<std::tuple<QString, QString, double>>> stnData;
    // 是否读到主测站行
    bool s = false;
    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine().trimmed();
        if (line == "") {
            continue;
        }

        QStringList lineList = line.split(",");
        if (lineList.size() != 3 && lineList.size() != 1) {
            qDebug() << "In2数据有误" << line;
            return std::nullopt;
        }

        if (!s) {
            if (lineList.size() == 1) {
                QString stn = lineList[0];
                stnData.first = stn;
                s = true;
            }
        } else {
            if (lineList.size() == 3) {
                stnData.second.append(std::make_tuple(lineList[0], lineList[1], lineList[2].toDouble()));
            } else if (lineList.size() == 1) {
                in.seek(pos);
                break;
            }
        }
    }
    return stnData;
}

QPair<QMap<QString, QPair<Eigen::Vector2d, bool>>, QList<QPair<QString, QString>>> QAdjustFileParse::Gra::ParseGra(QTextStream& in)
{
    QList<QPair<QString, QString>> lines;
    QMap<QString, QPair<Eigen::Vector2d, bool>> points;

    int pointNum = 0;
    int pointLineNum = 0;
    while (!in.atEnd()) {
        QString line = in.readLine();
        if (line.isEmpty()) {
            continue;
        }
        QStringList list = line.trimmed().split(",");
        if (list.size() == 4) {
            auto&& pname = list[0];
            auto&& x = list[1].toDouble();
            auto&& y = list[2].toDouble();
            bool isKnown = list[3].toInt();
            points[pname] = qMakePair(Eigen::Vector2d(x, y), isKnown);
        } else if (list.size() == 2) {
            lines.append(qMakePair(list[0], list[1]));
        }
    }
    return qMakePair(points, lines);
}

std::optional<QList<GSISeg>> QAdjustFileParse::GSI::ParseGSI2Entity(QTextStream& stream)
{
    QList<GSISeg> segs;
    auto&& opt = QAdjustFileParse::GSI::ParseGSI(stream);
    if (!opt.has_value()) {
        return std::nullopt;
    }
    auto&& gsis = opt.value();
    for (auto&& seg : gsis) {
        GSISeg gsiSeg;
        gsiSeg.setStartPoint(seg.second.first);
        gsiSeg.setAltitude(seg.second.second);

        QList<GSIRound> rounds;
        for (auto&& round : seg.first) {
            GSIRound gsiRound(round);
            rounds.append(gsiRound);
        }
        gsiSeg.setEndPoint(rounds.last().fName());
        gsiSeg.setGSIRounds(rounds);
        gsiSeg.setEndPoint(rounds.last().fName());
        gsiSeg.setDistance(rounds.last().distFromStart());
        gsiSeg.setHeightDiff(rounds.last().altitude() - gsiSeg.altitude());
        segs.append(gsiSeg);
    }
    return segs;
}
