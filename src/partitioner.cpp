#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include <utility>
#include <chrono>
#include <typeinfo>
#include <type_traits>
#include <algorithm>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


// Time measure utility source from: https://stackoverflow.com/questions/2808398/easily-measure-elapsed-time
template <
    class result_t   = chrono::milliseconds,
    class clock_t    = chrono::steady_clock,
    class duration_t = chrono::milliseconds
>
auto since(chrono::time_point<clock_t, duration_t> const& start) -> decltype( chrono::duration_cast<result_t>(clock_t::now() - start) )
{
    return chrono::duration_cast<result_t>(clock_t::now() - start);
}

int cmp(pair<int, int> a, pair<int, int> b) {
    return a.first < b.first;
}

void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }

    cout << endl;
    //cout << "==================== Parse Input ====================" << endl;
    cout << "[Info]: Load circuit with " << _cellNum << " cells and " << _netNum << " nets." << endl;
    cout << "        Balance factor: " << _bFactor << endl;

    return;
}

bool Partitioner::isBalanced(int moveCellId)
{
    int partAfterMove[2];
    partAfterMove[0] = _partSize[0] + (_cellArray[moveCellId]->getPart() ? 1 : -1); // if before move is B, than increase part A
    partAfterMove[1] = _partSize[1] + (_cellArray[moveCellId]->getPart() ? -1 : 1);

    // balance constraint
    if ( (1-_bFactor)/2 * _cellNum > partAfterMove[0] || (1+_bFactor)/2 * _cellNum < partAfterMove[0] ) {
        //cout << "[Info]: ... balance constraint doesn't satisfied with Part(A): "<< partAfterMove[0] <<" with range " << ((1-_bFactor)/2 * _cellNum) << " ~ " << ((1+_bFactor)/2 * _cellNum) << endl;
        return false;
    }
    if ( (1-_bFactor)/2 * _cellNum > partAfterMove[1] || (1+_bFactor)/2 * _cellNum < partAfterMove[1] ) {
        //cout << "[Info]: ... balance constraint doesn't satisfied with Part(B): "<< partAfterMove[1] <<" with range " << ((1-_bFactor)/2 * _cellNum) << " ~ " << ((1+_bFactor)/2 * _cellNum) << endl;
        return false;
    }
    return true;
}

int Partitioner::F(int netId, Cell* refCell)
{
    return _netArray[netId]->getPartCount( refCell->getPart() );
}

int Partitioner::T(int netId, Cell* refCell)
{
    return _netArray[netId]->getPartCount( !( refCell->getPart() ) );
}

void Partitioner::insertToBList(Cell* cell)
{
    // Auto detect which side by reference from `cell->getPart()`
    // Insert to the front of linked list
    Node* head = _bList[cell->getPart()][cell->getGain()];
    if (head != nullptr) {
        head->setPrev(cell->getNode());
        cell->getNode()->setNext(head);
    }
    _bList[cell->getPart()][cell->getGain()] = cell->getNode();
    return;
}

void Partitioner::removeFromBList(Cell* cell)
{
    // Remove cell from bucket list
    if (_bList[cell->getPart()][cell->getGain()] == cell->getNode()) { // if is the first node
        if (cell->getNode()->getNext() != nullptr)
            _bList[cell->getPart()][cell->getGain()] = cell->getNode()->getNext();
        else {
            _bList[cell->getPart()].erase( cell->getGain() ); // remove empty entry
            assert(cell != nullptr);
        }
    }
    // Remove from doubly linked list
    Node* prev = cell->getNode()->getPrev();
    if (prev != nullptr)
        prev->setNext( cell->getNode()->getNext() );
    Node* next = cell->getNode()->getNext();
    if (next != nullptr)
        next->setPrev(prev);

    cell->getNode()->setPrev(nullptr);
    cell->getNode()->setNext(nullptr);
    return;
}

void Partitioner::partition()
{
    const auto start = std::chrono::steady_clock::now();
    // Find max p(i)
    for (Cell* cell: _cellArray)
        if (_maxPinNum < cell->getPinNum())
            _maxPinNum = cell->getPinNum();

    // Init partition
    vector<pair<int, int>> pins;    // <pins, netId>
    for (Net* net: _netArray)
        pins.push_back( make_pair(net->getCellList().size(), _netName2Id[net->getName()]) );
    assert(pins.size() == _netArray.size());
    sort(pins.begin(), pins.end(), cmp);
    _partSize[0] = 0;
    _partSize[1] = 0;
    //srand(114514);
    //bool flip;
    for (auto x: pins) {
        for (int cellId: _netArray[x.second]->getCellList()) {
            if (_cellArray[cellId]->getLock())
                continue;

            if (_partSize[0] >= _cellNum * 0.5) {
                _cellArray[cellId]->move();
                _partSize[1]++;
            }
            else if (_partSize[1] >= _cellNum * 0.5) {
                _partSize[0]++;
            }
            else {
                _partSize[0]++;
            }

            _cellArray[cellId]->lock();
        }
    }
    /*
    for(int i = _cellNum * 0.5; i < _cellNum; i++)
        _cellArray[i]->move();
    */

    /* Update initial stats */
    // Count net->_partCount
    for (Net* net: _netArray) {
        net->setPartCount(0, 0);
        net->setPartCount(1, 0);
        for (int cellId: net->getCellList())
            net->incPartCount( _cellArray[cellId]->getPart() );
    }
    // Count _cutSize
    _initCutSize = 0;
    for (Net* net: _netArray)
        if (net->getPartCount(0) != 0 && net->getPartCount(1) != 0)
            _initCutSize++;
    //cout << "[Info]: Initial cutsize = " << _initCutSize << endl;

    // Fire optimization!!
    _iterNum = 0;
    do {
        //cout << endl;
        //cout << "==================== Iteration " << _iterNum <<" ====================" << endl;
        /* Update informations base on Cell states */
        // Reset everything except Cell->_part
        _partSize[0] = 0;
        _partSize[1] = 0;
        for (int i=0; i<2; i++)
            _bList[i].clear();
        // Init net->_partCount
        for (Net* net: _netArray) {
            net->setPartCount(0, 0);
            net->setPartCount(1, 0);
            for (int cellId: net->getCellList())
                net->incPartCount( _cellArray[cellId]->getPart() );
        }
        // Compute iteration initial gain, _partSize, _bList
        for (Cell* cell: _cellArray) {
            cell->setGain(0);
            cell->unlock();  // all cells are free at the begining of each iteration
            for (int netId: cell->getNetList()) {
                if (F(netId, cell) == 1) // F(n) == 1
                    cell->incGain();
                if (T(netId, cell) == 0) // T(n) == 0
                    cell->decGain();
            }
            // Insert to the front of linked list
            insertToBList(cell);
            // Update _partSize
            _partSize[cell->getPart()]++;
        }

        //cout << "[Info]: partSize(A) = " << _partSize[0] << ", partSize(B) = " << _partSize[1] << endl;

        _unlockNum[0] = _partSize[0];
        _unlockNum[1] = _partSize[1];
        vector<pair<int, int>> steps;   // (cellId, gain_head)
        while (_unlockNum[0] + _unlockNum[1] > 0) {
            // Find valid _maxGainCell
            _maxGainCell = nullptr;
            auto iter0 = _bList[0].rbegin();
            auto iter1 = _bList[1].rbegin();
            do {
                if ( iter0 == _bList[0].rend() ) {
                    _maxGainCell = iter1->second;
                    iter1++;
                } else if ( iter1 == _bList[1].rend() ) {
                    _maxGainCell = iter0->second;
                    iter0++;
                } else {
                    if ( iter0->first > iter1->first ) {
                        _maxGainCell = iter0->second;
                        iter0++;
                    } else {
                        _maxGainCell = iter1->second;
                        iter1++;
                    }
                }
                assert(_maxGainCell != nullptr);
            } while (!isBalanced(_maxGainCell->getId()) && ( iter0 != _bList[0].rend() || iter1 != _bList[1].rend() ) );
            assert(_maxGainCell != nullptr);
            assert(isBalanced(_maxGainCell->getId()));

            //cout << "[Info]: Select _maxGainCell " << _maxGainCell->getId() << " with gain "<< _cellArray[_maxGainCell->getId()]->getGain() << endl;

            Cell* selectedCell = _cellArray[_maxGainCell->getId()];
            // Lock cell
            selectedCell->lock();
            // Update gain
            for (int netId: selectedCell->getNetList()) {

                int cnt_t1 = 0, cnt_f1 = 0;
                if (T(netId, selectedCell) == 0) { // T(n) == 0 before move
                    for (int cellId: _netArray[netId]->getCellList())
                        if ( !(_cellArray[cellId]->getLock()) ) {
                            // remove from current bucket list
                            removeFromBList(_cellArray[cellId]);
                            // increase gain
                            _cellArray[cellId]->incGain();
                            // update position in bucket list
                            insertToBList(_cellArray[cellId]);
                        }
                }
                else if (T(netId, selectedCell) == 1) { // T(n) == 1 before move
                    for (int cellId: _netArray[netId]->getCellList())
                        if ( !(_cellArray[cellId]->getLock()) && (_cellArray[cellId] != selectedCell) && (_cellArray[cellId]->getPart() != selectedCell->getPart()) ) {
                            assert(cellId != _maxGainCell->getId());
                            assert(cnt_t1 < 1);
                            // remove from current bucket list
                            removeFromBList(_cellArray[cellId]);
                            // decrease gain
                            _cellArray[cellId]->decGain();
                            // update position in bucket list
                            insertToBList(_cellArray[cellId]);
                            cnt_t1++;
                            break;
                        }
                }
                // Update partCount to represent movement
                _netArray[netId]->decPartCount( selectedCell->getPart() );
                _netArray[netId]->incPartCount( !(selectedCell->getPart()) );

                if (F(netId, selectedCell) == 0) { // F(n) == 0 after move
                    for (int cellId: _netArray[netId]->getCellList())
                        if ( !(_cellArray[cellId]->getLock()) ) {
                            // remove from current bucket list
                            removeFromBList(_cellArray[cellId]);
                            // increase gain
                            _cellArray[cellId]->decGain();
                            // update position in bucket list
                            insertToBList(_cellArray[cellId]);
                        }
                }
                else if (F(netId, selectedCell) == 1) { // F(n) == 1 after move
                    for (int cellId: _netArray[netId]->getCellList())
                        if ( !(_cellArray[cellId]->getLock()) && (_cellArray[cellId] != selectedCell) && (_cellArray[cellId]->getPart() == selectedCell->getPart()) ) {
                            assert(cellId != _maxGainCell->getId());
                            assert(cnt_f1 < 1);
                            // remove from current bucket list
                            removeFromBList(_cellArray[cellId]);
                            // decrease gain
                            _cellArray[cellId]->incGain();
                            // update position in bucket list
                            insertToBList(_cellArray[cellId]);
                            cnt_f1++;
                            break;
                        }
                }
            }

            // Update _unlockNum
            _unlockNum[selectedCell->getPart()]--;
            // Remove cell from bucket list
            removeFromBList(selectedCell);
            //cout << "[Info]: ... lock _maxGainCell with gain "<< selectedCell->getGain() << " remain " << (_unlockNum[0] + _unlockNum[1]) << " unlocked cell waiting for process." << endl;
            // Update _partSize to refer the movement
            _partSize[selectedCell->getPart()]--;
            _partSize[!(selectedCell->getPart())]++;
            // Store move
            steps.push_back( make_pair(_maxGainCell->getId(), selectedCell->getGain()) );
        }

        // Find _maxAccGain
        _accGain = 0;
        _moveNum = -1;
        _maxAccGain = 0;
        _bestMoveNum = -1;
        for (auto const& step: steps) {
            _accGain += step.second;
            _moveNum++;

            if (_bestMoveNum < 0 || _maxAccGain < _accGain) {
                /* TODO: Choose move best to balance */

                _maxAccGain = _accGain;
                _bestMoveNum = _moveNum;
            }
        }
        // Move cells
        if (_maxAccGain > 0) {
            _moveNum = 0;
            for (auto const& step: steps) {
                if (_moveNum > _bestMoveNum)
                    break;
                _cellArray[step.first]->move();
                _moveNum++;
            }
            /*
            // Count net->_partCount
            for (Net* net: _netArray) {
                net->setPartCount(0, 0);
                net->setPartCount(1, 0);
                for (int cellId: net->getCellList())
                    net->incPartCount( _cellArray[cellId]->getPart() );
            }
            // Count _cutSize
            _cutSize = 0;
            for (Net* net: _netArray)
                if (net->getPartCount(0) != 0 && net->getPartCount(1) != 0)
                    _cutSize++;
            cout << "[Info]: Move " << _bestMoveNum + 1 << " cell(s)." << endl;
            cout << "[Info]: Optmized cutsize: " << _cutSize << endl;
            cout << "[Info]: Max acc-gain: " << _maxAccGain << endl;
            */
        }

        _iterNum++;

    } while (_maxAccGain > 0);

    /* Update partition results */
    // Count net->_partCount
    for (Net* net: _netArray) {
        net->setPartCount(0, 0);
        net->setPartCount(1, 0);
        for (int cellId: net->getCellList())
            net->incPartCount( _cellArray[cellId]->getPart() );
    }
    // Count _partSize
    _partSize[0] = 0;
    _partSize[1] = 0;
    for (Cell* cell: _cellArray)
        _partSize[cell->getPart()]++;
    // Count _cutSize
    _cutSize = 0;
    for (Net* net: _netArray)
        if (net->getPartCount(0) != 0 && net->getPartCount(1) != 0)
            _cutSize++;

    cout << "[Info]: Elapsed (ms): " << since(start).count() << std::endl;

    return;
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Intial cutsize:    " << _initCutSize << endl;
    cout << " Optmized Cutsize:  " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
