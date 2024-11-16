#include "QAdjustFileParse.h"
#include "QAdjustFileParseImpl.h"

#include <QDebug>
#include <QFile>
#include <QJsonArray>
#include <QMap>
#include <QPair>
#include <QTextStream>

#include <regex>

#include "DatRound.h"
#include "GSIRound.h"
#include "QtUtils.h"

using namespace QAdjustFileParse;

QList<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>> QAdjustFileParse::Dat::ParseDat(QTextStream& in)
{
    QList<QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>>> dat;
    while (!in.atEnd()) {
        dat.push_back(Dat::ParseDatSeg(in));
    }
    return dat;
}

QADJUSTFILEPARSE_EXPORT QList<DatSeg> QAdjustFileParse::Dat::ParseDat2Entity(QTextStream& stream)
{
    QList<DatSeg> segs;
    auto&& dat = QAdjustFileParse::Dat::ParseDat(stream);
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

QPair<QList<QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>>>, QPair<QString, double>> QAdjustFileParse::Dat::ParseDatSeg(QTextStream& in)
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
        if (line.contains("Start-Line")) {
            s = true;
            start = Dat::ParseDatRoundAltitudeLine(in.readLine());
        } else if (line.contains("End-Line")) {
            s = false;
            break;
        } else {
            in.seek(pos);
        }

        if (s) {
            auto&& data = Dat::ParseDatRound(in);
            seg.push_back(data);
        }
    }
    return qMakePair(seg, start);
}

QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>> QAdjustFileParse::Dat::ParseDatRound(QTextStream& in)
{
    QPair<QList<std::tuple<QString, QString, double, double>>, QPair<QString, double>> round;
    std::optional<QString> startPoint { std::nullopt };
    while (!in.atEnd()) {
        QString line = in.readLine();
        if (!Dat::IsValidDatLine(line)) {
            continue;
        }
        auto&& array = Dat::ParseDatRoundLine(line);
        if (!startPoint.has_value()) {
            startPoint = std::get<0>(array);
            round.first.push_back(array);
        } else if (startPoint == std::get<0>(array)) {
            round.first.push_back(array);
            round.second = ParseDatRoundAltitudeLine(in.readLine());
            break;
        } else {
            round.first.push_back(array);
        }
    }
    return round;
}

std::tuple<QString, QString, double, double> QAdjustFileParse::Dat::ParseDatRoundLine(QString line)
{
    QStringList parts = line.split("|");
    parts.removeFirst();
    parts.removeFirst();
    QString target, orient;
    double heightDiff, distance;
    target = SpaceReplacSplit(parts[0])[1];
    orient = SpaceReplacSplit(parts[1])[0];
    heightDiff = SpaceReplacSplit(parts[1])[1].toDouble();
    distance = SpaceReplacSplit(parts[2])[1].toDouble();
    return std::make_tuple(target, orient, heightDiff, distance);
}

QPair<QString, double> QAdjustFileParse::Dat::ParseDatRoundAltitudeLine(QString line)
{
    QStringList parts = line.split("|");
    parts.removeFirst();
    parts.removeFirst();
    QString target = SpaceReplacSplit(parts[0])[1];
    double altitude = SpaceReplacSplit(parts[3])[1].toDouble();
    return qMakePair(target, altitude);
}

QMap<QPair<QString, QString>, QPair<double, double>> QAdjustFileParse::In1::ParseIn1(QTextStream& in)
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
        if (line.split(",").size() == 4) {
            s = true;
        }
        if (!s) {
            // QStringList splitString = line.split(",");
            // QString pname = splitString[0];
            // double altitude = splitString[1].toDouble();
            // knownPoints[pname] = altitude;
        } else {
            QStringList splitString = line.split(",");
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
        }
    }

    return data;
}

QList<In1Observed> QAdjustFileParse::In1::ParseIn12Enity(QTextStream& stream)
{
    QList<In1Observed> obss;
    auto&& in1 = QAdjustFileParse::In1::ParseIn1(stream);

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
    if (line.contains("####") || line.contains("Measurement repeated") || line.contains("Station repeated")) {
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

QList<QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>>> QAdjustFileParse::GSI::ParseGSI(QTextStream& in)
{
    QList<QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>>> gsi;
    while (!in.atEnd()) {
        gsi.push_back(GSI::ParseGSISeg(in));
    }
    return gsi;
}

QADJUSTFILEPARSE_EXPORT QList<GSISeg> QAdjustFileParse::GSI::ParseGSI2Entity(QTextStream& stream)
{
    QList<GSISeg> segs;
    auto&& gsis = QAdjustFileParse::GSI::ParseGSI(stream);

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

QPair<QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>>, QPair<QString, double>> QAdjustFileParse::GSI::ParseGSISeg(QTextStream& in)
{
    QPair<QString, double> start;
    QList<QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>>> rounds;
    bool s = false;
    while (!in.atEnd()) {
        auto&& pos = in.pos();
        QString line = in.readLine();
        QStringList list = SpaceReplacSplit(line);
        if (s) {
            in.seek(pos);
            if (list.size() == 1) {
                break;
            } else {
                rounds.append(GSI::ParseGSIRound(in));
                continue;
            }
        }

        if (list.size() == 1) {
            s = true;
            // 处理起始点数据
            QString startLine = in.readLine();
            QStringList startList = SpaceReplacSplit(startLine);
            start.first = GSI::ParseGSIPName(startList[0]);
            start.second = GSI::ParseGSIValue(startList[1]);
            continue;
        }
    }
    return qMakePair(rounds, start);
}

QPair<QList<std::tuple<QString, bool, double, double>>, std::tuple<QString, double, double, double, double, double>> QAdjustFileParse::GSI::ParseGSIRound(QTextStream& in)
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

    for (auto&& rLine : rLines) {
        QString target = GSI::ParseGSIPName(rLine[0]);
        double rulerNum = GSI::ParseGSIValue(rLine[1]);
        double dist = GSI::ParseGSIValue(rLine[2]);
        // 前进目标方向为true
        bool orient = false;
        if (target == targetPName) {
            orient = true;
        }
        rdatas.append(std::make_tuple(target, orient, rulerNum, dist));
    }
    return qMakePair(rdatas, sdata);
}

QString QAdjustFileParse::GSI::ParseGSIPName(QString input)
{
    return input.right(8);
}

double QAdjustFileParse::GSI::ParseGSIValue(QString input)
{
    double data = input.right(8).toInt() / 100000.0;
    if (input.contains("-")) {
        return -data;
    }
    return data;
}

QPair<QString, QList<QPair<int, QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>>>>> QAdjustFileParse::SUC::ParseSuc(QTextStream& in)
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
        }
        if (s) {
            int surveyTime = lineList[0].trimmed().toInt();
            rounds.append(qMakePair(surveyTime, SUC::ParseSucRound(in)));
        }
    }
    return qMakePair(deviceName, rounds);
}

QList<SucPoint> QAdjustFileParse::SUC::ParseSuc2Entity(QTextStream& stream)
{
    QList<SucPoint> points;
    auto&& suc = QAdjustFileParse::SUC::ParseSuc(stream);
    QString stnName = suc.first;
    for (auto&& round : suc.second) {
        int surveyTime = round.first;
        for (auto&& point : round.second) {
            QString pname = point.first;
            auto&& [hf, vf, sf] = std::get<0>(point.second).first;
            auto&& [hr, vr, sr] = std::get<0>(point.second).second;
            auto&& deviceHeight = std::get<1>(point.second);
            auto&& prismHeight = std::get<2>(point.second);
            SucPoint point(stnName, surveyTime, pname, Utils::HV::AMS2A(hf), Utils::HV::AMS2A(vf), sf, Utils::HV::AMS2A(hr), Utils::HV::AMS2A(vr), sr, deviceHeight, prismHeight);
            points.append(point);
        }
    }
    return points;
}

QList<QPair<QString, std::tuple<QPair<std::tuple<double, double, double>, std::tuple<double, double, double>>, double, double>>> QAdjustFileParse::SUC::ParseSucRound(QTextStream& in)
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
        }
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

QList<QPair<QString, QList<std::tuple<QString, QString, double>>>> QAdjustFileParse::In2::ParseIn2(QTextStream& in)
{
    // 已知误差
    // std::tuple<double, double, double, double> knownDeviation;
    // 已知测站数据 {测站名，{x,y}}
    // QMap<QString, Eigen::Vector2d> knownStns;
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
            // 是否是误差行
            // if (Utils::IsNumber(lineList[0])) {
            //    if (lineList.size() == 3) {
            //        //knownDeviation = std::make_tuple(lineList[0].toDouble(), 0, lineList[1].toDouble(), lineList[2].toDouble());
            //    } else if (lineList.size() == 4) {
            //        knownDeviation = std::make_tuple(lineList[0].toDouble(), lineList[1].toDouble(), lineList[2].toDouble(), lineList[3].toDouble());
            //    }
            //} else {
            //    if (lineList.size() == 3) {
            //        knownStns[lineList[0]] = Eigen::Vector2d(lineList[1].toDouble(), lineList[2].toDouble());
            //    }
            //}
        }

        if (lineList.size() == 1) {
            in.seek(pos);
            stnDatas.append(ParseIn2StnData(in));
            s = true;
        }
    }
    return stnDatas;
}

QList<QPair<QString, QList<In2Observed>>> QAdjustFileParse::In2::ParseIn22Entity(QTextStream& stream)
{
    auto&& ObservedGroups = QAdjustFileParse::In2::ParseIn2(stream);
    QList<QPair<QString, QList<In2Observed>>> result;
    for (auto&& ObservedGroup : ObservedGroups) {
        QPair<QString, QList<In2Observed>> stnData;
        QList<In2Observed> obss;
        QString stn = ObservedGroup.first;
        for (auto&& obsData : ObservedGroup.second) {
            auto&& [pname, type, value] = obsData;
            In2Observed obs(stn, pname, S2E(In2Observed::ObservedValueType, type), value);
            obss.append(obs);
        }
        stnData = qMakePair(stn, obss);
        result.append(stnData);
    }
    return result;
}

QPair<QString, QList<std::tuple<QString, QString, double>>> QAdjustFileParse::In2::ParseIn2StnData(QTextStream& in)
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

std::tuple<QMap<QPair<QString, QString>, bool>, QMap<QString, Eigen::Vector2d>, QMap<QString, std::tuple<double, double, double>>> QAdjustFileParse::Gra::ParseGra(QTextStream& in)
{
    QMap<QPair<QString, QString>, bool> lines;
    QMap<QString, Eigen::Vector2d> points;
    QMap<QString, std::tuple<double, double, double>> ovals;
    int pointNum = 0;
    int pointLineNum = 0;
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList list = SpaceReplacSplit(line);
        if (list.size() == 2) {
            pointNum = list[0].toInt();
            continue;
        }

        if (pointNum != 0) {
            // 读取点行数据
            if (pointLineNum < pointNum) {
                if (list.size() == 3) {
                    QString pname = list[0];
                    double x = list[1].toDouble();
                    double y = list[2].toDouble();
                    points[pname] = Eigen::Vector2d(x, y);
                    pointLineNum++;
                } else if (list.size() == 6) {
                    QString pname = list[0];
                    double x = list[1].toDouble();
                    double y = list[2].toDouble();
                    // 长半轴半径
                    double a = list[3].toDouble();
                    // 短半轴半径
                    double b = list[4].toDouble();
                    // 长半轴旋转角
                    double a_angle = list[5].toDouble();
                    points[pname] = Eigen::Vector2d(x, y);
                    ovals[pname] = std::make_tuple(a, b, a_angle);
                    pointLineNum++;
                }
                continue;
            }
            // 读取行数据
            if (list.size() == 3) {
                QString p1 = list[0];
                QString p2 = list[1];
                bool isKnownLine = list[2] == "S" ? true : false;
                lines[qMakePair(p1, p2)] = isKnownLine;
            }
        }
    }
    return std::make_tuple(lines, points, ovals);
}
