/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_MAPS_H
#define MECHSYS_MAPS_H

// Std lib
#include <cstring>  // for strcmp
#include <iostream> // for cout
#include <sstream>  // for istringstream, ostringstream
#include <cstdarg>  // for va_list, va_start, va_end
#include <cmath>    // for fabs
#include <fstream>
#include <map>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/string.h>
#include <mechsys/util/numstreams.h>


/////////////////////////////////////////////////////////////////////////////////////////// SDPair


typedef std::map<String,double> StrDbl_t;

/** String-Double Pair. */
class SDPair : public StrDbl_t
{
public:
    // Constructors
    SDPair () {}                               ///< Default constructor
    SDPair (SDPair const & R) { (*this) = R; } ///< Copy constructor

    // Set methods
    void Set (const char * Str, ...); ///< Set given keys and vals: Set("ux uy", 1.0,2.0). NOTE: vals must be double
    void SetZero (Array<String> const & keys); ///< Create pair with all keys set to zero

    // Operators
    double       & operator() (char   const * Key);
    double const & operator() (char   const * Key) const;
    double       & operator() (String const & Key)       { return operator()(Key.CStr()); }
    double const & operator() (String const & Key) const { return operator()(Key.CStr()); }
    void           operator=  (SDPair const & R); ///< Assignment operator
    void           operator+= (SDPair const & R); ///< Used to merge keys, but not sum them up

    // Methods
    void   Del       (const char * Key);             ///< Delete key (and value) from pair
    size_t ReSet     (const char * Key, double Val); ///< Re-set value (creates if it doesn't exist). Return position in Keys of changed/set value 
    size_t AddVal    (const char * Key, double Val); ///< Add value to existent one (creates if it doesn't exist). Return position in Keys of changed/set value 
    void   SetValues (double Val);                   ///< Set all values equal to Val
    double ValOrZero (char const   * Key) const;     ///< Returns value for Key if Key exists, otherwise returns zero
    double ValOrZero (String const & Key) const { return ValOrZero(Key.CStr()); }
    bool   HasKey    (char const   * Key) const;
    bool   HasKeyLen (char const   * Key, int KeyLen) const;
    void   GetSubKeys(char const   * Key, int KeyLen, Array<String> & SubKeys) const;
    bool   HasKey    (String const & Key) const;
    long   IdxKey    (String const & Key) const { return Keys.Find(Key); } ///< Returns the index of key if found, otherwise returns -1
    void   Val2Key   (double Val, String & Key, double Tol=1.0e-15) const;
    void   clear     () { Keys.Clear(); StrDbl_t::clear(); }

    // Data
    Array<String> Keys;

#ifdef USE_BOOST_PYTHON
    void PySet (BPy::dict const & Pairs);
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// SIPair


typedef std::map<String,int> StrInt_t;

/** String-Int Pair. */
class SIPair : public StrInt_t
{
public:
    // Set methods
    void Set (const char * Str, ...); ///< Set given keys and vals: Set("a b", 1,2). NOTE: vals must be int

    // Operators
    int       & operator() (char const   * Key);
    int const & operator() (char const   * Key) const;
    int       & operator() (String const & Key)        { return operator()(Key.CStr()); }
    int const & operator() (String const & Key) const  { return operator()(Key.CStr()); }

    // Methods
    void Del    (const char * Key);             ///< Delete key (and value) from pair
    bool HasKey (char const   * Key) const;
    bool HasKey (String const & Key) const;
    void clear  () { Keys.Clear(); StrInt_t::clear(); }

    // Data
    Array<String> Keys;
};


/////////////////////////////////////////////////////////////////////////////////////////// Dict


typedef std::map<int, SDPair> Dict_t;

/** Dictionary. */
class Dict : public Dict_t
{
public:
    // Constructors
    Dict () {}                             ///< Default constructor
    Dict (Dict const & R) { (*this) = R; } ///< Copy constructor

    // Methods
    void Set     (int Key, const char * Str, ...);      ///< Set given key, keys, and vals: Set(-1, "ux uy", 1.0,2.0). NOTE: valus must be double
    void Set     (int Key, SDPair const & P);           ///< Set given SDPair (does merge)
    void SetZero (int Key, Array<String> const & keys); ///< Create dictionary with all keys set to zero
    void Del     (int Key);                             ///< Delete item corresponding to Key

    // Operators
    SDPair const & operator() (int Key) const;
    SDPair       & operator() (int Key);
    void           operator=  (Dict const & R); ///< Assignment operator
    void           operator+= (Dict const & R); ///< Used to merge keys

    // Methods
    bool HasKey (int Key) const;
    void clear  () { Keys.Clear(); Dict_t::clear(); }

    // Data
    Array<int> Keys;

#ifdef USE_BOOST_PYTHON
    void PySet (int Key, BPy::dict const & Pairs);
    void PyGet (int Key, BPy::dict       & Pairs);
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Table


typedef std::map<String,Array<double> > Table_t;

class Table : public Table_t
{
public:
    // Constructor
    Table () : NRows(0) {}

    // Methods
    /** Ex: Set("ux uy", 2, 1.0, 2.0,
     *                      3.0, 4.0);  */
    void Set     (const char * StrKeys, size_t NumRows, ...);
    void SetZero (const char * StrKeys, size_t NumRows);
    void Read    (const char * FileName);
    void Write   (const char * FileName, char const * NF="%15.8e");

    // Operators
    Array<double>       & operator() (String const & Key);
    Array<double> const & operator() (String const & Key) const;
    double              & operator() (String const & Key, size_t iRow);
    double        const & operator() (String const & Key, size_t iRow) const;

    // Data
    size_t        NRows;
    Array<String> Keys;

#ifdef USE_BOOST_PYTHON
    void PySet (BPy::list const & StrKeys, BPy::list const & Values);
#endif
};


///////////////////////////////////////////////////////////////////////////////////// operator << //////////////


std::ostream & operator<< (std::ostream & os, SDPair const & P)
{
    int nkeys = P.size();
    int k     = 0;
    /*
    os << "{";
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        String key = P.Keys[i];
        String val;  val.Printf("%g",P(key));
        os << "'" << key << "':" << val;
        if (k<nkeys-1) os << ", ";
        k++;
    }
    os << "}";
    */
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        String key = P.Keys[i];
        String val;  val.Printf("%g",P(key));
        os << key << "=" << val;
        if (k<nkeys-1) os << "  ";
        k++;
    }
    return os;
}

std::ostream & operator<< (std::ostream & os, SIPair const & P)
{
    int nkeys = P.size();
    int k     = 0;
    os << "{";
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        String key = P.Keys[i];
        String val;  val.Printf("%d",P(key));
        os << "'" << key << "':" << val;
        if (k<nkeys-1) os << ", ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, Dict const & D)
{
    int nkeys = D.size();
    int k = 0;
    os << "{";
    for (size_t i=0; i<D.Keys.Size(); ++i)
    {
        int key = D.Keys[i];
        os << key << ":" << D(key);
        if (k<nkeys-1) os << ",\n ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, Table const & T)
{
    // keys
    for (size_t i=0; i<T.Keys.Size(); ++i) os << Util::_8s << T.Keys[i];
    os << "\n";

    // values
    for (size_t i=0; i<T.NRows; ++i)
    {
        for (size_t j=0; j<T.Keys.Size(); ++j) os << Util::_8s << T(T.Keys[j],i);
        os << "\n";
    }

    return os;
}


/////////////////////////////////////////////////////////////////////////////////// SDPair: Implementation /////


inline void SDPair::Set(const char * Str, ...)
{
    std::istringstream iss(Str);
    String    key;
    va_list   arg_list;
    va_start (arg_list, Str);
    while (iss>>key)
    {
        Keys.XPush (key);
        (*this)[key] = va_arg(arg_list,double);
    }
    va_end (arg_list);
}

inline void SDPair::SetZero (Array<String> const & keys)
{
    for (size_t i=0; i<keys.Size(); ++i)
    {
        Keys.XPush (keys[i]);
        (*this)[keys[i]] = 0.0;
    }
}

inline double & SDPair::operator() (char const * Key)
{
    StrDbl_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::operator(): String-Double pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    return p->second;
}

inline double const & SDPair::operator() (char const * Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::operator(): String-Double pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    return p->second;
}

inline void SDPair::operator= (SDPair const & R)
{
    this->clear();
    Keys.Resize (R.Keys.Size());
    for (size_t i=0; i<R.Keys.Size(); ++i)
    {
        Keys[i]          = R.Keys[i];
        (*this)[Keys[i]] = R(Keys[i]);
    }
}

inline void SDPair::operator+= (SDPair const & R)
{
    for (size_t i=0; i<R.Keys.Size(); ++i)
    {
        if (this->HasKey(R.Keys[i]))
        {
            std::ostringstream oss;
            oss << (*this);
            throw new Fatal("SDPair::operator+= : SDPair [%s]\n   already contain key %s (values cannot be summed up).",oss.str().c_str(),R.Keys[i].CStr());
        }
        // This is fine though (if addition is permitted): (*this)(R.Keys[i]) += R(R.Keys[i]);
        else  this->Set(R.Keys[i].CStr(), R(R.Keys[i]));
    }
}

inline void SDPair::Del (const char * Key)
{
    StrDbl_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::Del: String-Double pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    else
    {
        Keys.DelVal (Key);
        this->erase (p);
    }
}

inline size_t SDPair::ReSet (const char * Key, double Val)
{
    StrDbl_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        (*this)[Key] = Val;
        Keys.Push (Key);
        return Keys.Size()-1;
    }
    else
    {
        p->second = Val;
        long pos = Keys.Find(Key);
        return static_cast<size_t>(pos);
    }
}

inline size_t SDPair::AddVal (const char * Key, double Val)
{
    StrDbl_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        (*this)[Key] = Val;
        Keys.Push (Key);
        return Keys.Size()-1;
    }
    else
    {
        p->second += Val;
        long pos = Keys.Find(Key);
        return static_cast<size_t>(pos);
    }
}

inline void SDPair::SetValues (double Val)
{
    for (StrDbl_t::iterator p=this->begin(); p!=this->end(); ++p)
        p->second = Val;
}

inline double SDPair::ValOrZero (char const * Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    if (p==this->end()) return 0.0;
    else                return p->second;
}

inline bool SDPair::HasKey (char const * Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}

inline bool SDPair::HasKeyLen (char const * Key, int KeyLen) const
{
    for (StrDbl_t::const_iterator p=this->begin(); p!=this->end(); ++p)
    {
        if (p->first.substr(0, KeyLen) == Key)
        {
            return true;
        }
    }
    return false;
}

// example: "mainkey_subkey1_subkey2" => SubKeys = [subkey1, subkey2]
inline void SDPair::GetSubKeys (char const * MainKey, int MainKeyLen, Array<String> & SubKeys) const
{
    for (StrDbl_t::const_iterator p=this->begin(); p!=this->end(); ++p)
    {
        if (p->first.substr(0, MainKeyLen) == MainKey)
        {
            std::stringstream subkeys(p->first.substr(MainKeyLen));
            std::string item;
            while (std::getline(subkeys, item, '_'))
            {
                if (item.size() > 0)
                {
                    SubKeys.Push(item);
                }
            }
            return;
        }
    }
}

inline bool SDPair::HasKey (String const & Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}

inline void SDPair::Val2Key (double Val, String & Key, double Tol) const
{
    bool found = false;
    for (StrDbl_t::const_iterator p=this->begin(); p!=this->end(); ++p)
    {
        if (fabs(Val-p->second)<Tol)
        {
            Key   = p->first;
            found = true;
            break;
        }
    }
    if (!found)
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::Val2Key: Could not find Val=%g in map: %s",Val,oss.str().c_str());
    }
}

#ifdef USE_BOOST_PYTHON
inline void SDPair::PySet (BPy::dict const & Pairs)
{
    BPy::object const & keys = Pairs.iterkeys();
    BPy::object const & vals = Pairs.itervalues();
    for (int i=0; i<BPy::len(Pairs); ++i)
    {
        char const * key = BPy::extract<char const *>(keys.attr("next")());
        double       val = BPy::extract<double      >(vals.attr("next")());
        Set (key, val);
    }
}
#endif


/////////////////////////////////////////////////////////////////////////////////// SIPair: Implementation /////


inline void SIPair::Set(const char * Str, ...)
{
    std::istringstream iss(Str);
    String    key;
    va_list   arg_list;
    va_start (arg_list, Str);
    while (iss>>key)
    {
        Keys.XPush (key);
        (*this)[key] = va_arg(arg_list,int);
    }
    va_end (arg_list);
}

inline int & SIPair::operator() (char const * Key)
{
    StrInt_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SIPair::operator(): String-Int pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    return p->second;
}

inline int const & SIPair::operator() (char const * Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SIPair::operator(): String-Int pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    return p->second;
}

inline void SIPair::Del (const char * Key)
{
    StrInt_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SIPair::Del: String-Int pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    else
    {
        Keys.DelVal (Key);
        this->erase (p);
    }
}

inline bool SIPair::HasKey (char const * Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}

inline bool SIPair::HasKey (String const & Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}


///////////////////////////////////////////////////////////////////////////////////// Dict: Implementation /////


inline void Dict::Set (int Key, const char * Str, ...)
{
    bool has_key = HasKey(Key);
    SDPair & pair = (has_key ? (*this)(Key) : (*this)[Key]);
    std::istringstream iss(Str);
    String    skey;
    va_list   arg_list;
    va_start (arg_list, Str);
    while (iss>>skey) pair.Set (skey.CStr(), va_arg(arg_list,double));
    va_end (arg_list);
    if (!has_key) Keys.Push (Key);
}

inline void Dict::Set (int Key, SDPair const & P)
{
    for (size_t i=0; i<P.Keys.Size(); ++i)
        this->Set (Key, P.Keys[i].CStr(), P(P.Keys[i]));
    if (P.Keys.Size()==0)
    {
        if (!HasKey(Key))
        {
            (*this)[Key] = P;
            Keys.Push (Key);
        }
    }
}

inline void Dict::SetZero (int Key, Array<String> const & keys)
{
    SDPair pair;
    pair.SetZero (keys);
    this->Set (Key, pair);
}

inline void Dict::Del (int Key)
{
    Dict_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("Dict::Del: Dictionary: %s does not have a key = %d",oss.str().c_str(),Key);
    }
    erase       (p);
    Keys.DelVal (Key);
}

inline SDPair const & Dict::operator() (int Key) const
{
    Dict_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("Dict::operator(): Dictionary: %s does not have a key = %d",oss.str().c_str(),Key);
    }
    return p->second;
}

inline SDPair & Dict::operator() (int Key)
{
    Dict_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("Dict::operator(): Dictionary: %s does not have a key = %d",oss.str().c_str(),Key);
    }
    return p->second;
}

inline void Dict::operator= (Dict const & R)
{
    this->clear();
    Keys.Resize (R.Keys.Size());
    for (size_t i=0; i<R.Keys.Size(); ++i)
    {
        Keys[i]          = R.Keys[i];
        (*this)[Keys[i]] = R(Keys[i]);
    }
}

inline void Dict::operator+= (Dict const & R)
{
    for (size_t i=0; i<R.Keys.Size(); ++i)
    {
        if (this->HasKey(R.Keys[i])) (*this)(R.Keys[i]) += R(R.Keys[i]); // merge
        else  this->Set(R.Keys[i], R(R.Keys[i]));
    }
}

inline bool Dict::HasKey (int Key) const
{
    if (Keys.Find(Key)<0) return false;
    else                  return true;
}

#ifdef USE_BOOST_PYTHON
inline void Dict::PySet (int Key, BPy::dict const & Pairs)
{
    BPy::object const & keys = Pairs.iterkeys();
    BPy::object const & vals = Pairs.itervalues();
    for (int i=0; i<BPy::len(Pairs); ++i)
    {
        char const * key = BPy::extract<char const *>(keys.attr("next")());
        double       val = BPy::extract<double      >(vals.attr("next")());
        Set (Key, key, val);
    }
}

inline void Dict::PyGet (int Key, BPy::dict & Pairs)
{
    SDPair const & pairs = (*this)(Key);
    for (size_t i=0; i<pairs.Keys.Size(); ++i) Pairs[pairs.Keys[i].CStr()] = pairs(pairs.Keys[i]);
}
#endif


//////////////////////////////////////////////////////////////////////////////////// Table: Implementation /////


inline void Table::Set (const char * StrKeys, size_t NumRows, ...)
{
    NRows = NumRows;
    Keys.Resize(0);

    // retrieve keys and initialize table
    String             key;
    std::istringstream iss(StrKeys);
    while (iss>>key)
    {
        if (Keys.Find(key)<0) Keys.Push(key);
        (*this)[key].Resize (NRows);
    }

    // read values
    va_list   arg_list;
    va_start (arg_list, NumRows);
    for (size_t i=0; i<NumRows; ++i)
    {
        for (size_t j=0; j<Keys.Size(); ++j)
        {
            (*this)[Keys[j]][i] = va_arg(arg_list,double);
        }
    }
    va_end (arg_list);
}

inline void Table::SetZero (const char * StrKeys, size_t NumRows)
{
    NRows = NumRows;
    Keys.Resize(0);

    // retrieve keys and initialize table
    String             key;
    std::istringstream iss(StrKeys);
    while (iss>>key)
    {
        if (Keys.Find(key)<0) Keys.Push(key);
        (*this)[key].Resize    (NRows);
        (*this)[key].SetValues (0.0);
    }
}

inline void Table::Read (const char * FileName)
{
    // open file
    std::fstream fil(FileName, std::ios::in);
    if (!fil.is_open()) throw new Fatal("Table::Read Could not open file < %s >",FileName);

    // erase data
    Table_t::clear();
    NRows = 0;

    // parse
    bool header = true;
    String line, str;
    double val;
    while (!fil.eof())
    {
        // read line
        std::getline (fil,line);
        std::istringstream iss(line);

        // header
        if (header)
        {
            while (iss>>str)
            {
                Keys.Push (str);
                (*this)[str].Resize (0);
            }
            header = false;
        }
        else
        {
            for (size_t i=0; i<Keys.Size(); ++i)
            {
                if (iss>>val)
                {
                    (*this)[Keys[i]].Push (val);
                    if (i==0) NRows++;
                }
            }
        }
    }
}

inline void Table::Write (const char * FileName, char const * NF)
{
    // keys
    String fmt, buf;
    fmt.TextFmt (NF);
    std::ostringstream oss;
    for (size_t i=0; i<Keys.Size(); ++i)
    {
        buf.Printf (fmt.CStr(), Keys[i].CStr());
        oss << buf << " ";
    }
    oss << "\n";

    // values
    for (size_t i=0; i<NRows; ++i)
    {
        for (size_t j=0; j<Keys.Size(); ++j)
        {
            buf.Printf (NF, operator()(Keys[j],i));
            oss << buf << " ";
        }
        oss << "\n";
    }

    // open file and save data
    std::ofstream of(FileName, std::ios::out);
    of << oss.str();
    of.close();
}

inline Array<double> & Table::operator() (String const & Key)
{
    Table_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << Keys;
        throw new Fatal("Table::operator(): Table with Keys=[%s] does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second;
}

inline Array<double> const & Table::operator() (String const & Key) const
{
    Table_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << Keys;
        throw new Fatal("Table::operator(): Table with Keys=[%s] does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second;
}

inline double & Table::operator() (String const & Key, size_t iRow)
{
    Table_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << Keys;
        throw new Fatal("Table::operator(): Table with Keys=[%s] does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second[iRow];
}

inline double const & Table::operator() (String const & Key, size_t iRow) const
{
    Table_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << Keys;
        throw new Fatal("Table::operator(): Table with Keys=[%s] does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second[iRow];
}

#ifdef USE_BOOST_PYTHON
inline void Table::PySet (BPy::list const & StrKeys, BPy::list const & Values)
{
    NRows = BPy::len(Values);

    // retrieve keys and initialize table
    size_t nkeys = BPy::len(StrKeys);
    for (size_t i=0; i<nkeys; ++i)
    {
        char const * key = BPy::extract<char const *>(StrKeys[i])();
        if (Keys.Find(key)<0) Keys.Push(key);
        (*this)[key].Resize (NRows);
    }

    // read values
    for (size_t i=0; i<NRows; ++i)
    {
        BPy::list const & line = BPy::extract<BPy::list>(Values[i])();
        for (size_t j=0; j<Keys.Size(); ++j)
        {
            (*this)[Keys[j]][i] = BPy::extract<double>(line[j])();
        }
    }
}
#endif


#endif // MECHSYS_MAPS_H
