// Copyright (c) 1996,1998 The Regents of the University of California.
// All rights reserved.
// 
// Redistribution and use in source and binary forms are permitted
// provided that the above copyright notice and this paragraph are
// duplicated in all such forms and that any documentation,
// distribution and/or use acknowledge that the software was developed
// by the Computer Graphics Laboratory, University of California,
// San Francisco.  The name of the University may not be used to
// endorse or promote products derived from this software without
// specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// IN NO EVENT SHALL THE REGENTS OF THE UNIVERSITY OF CALIFORNIA BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
// OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS SOFTWARE.

// $Id: Symbol.h 36238 2012-04-26 00:05:10Z goddard $

#ifndef chimera_Symbol_h
# define chimera_Symbol_h

# include "_chimera_config.h"
# include <iostream>
# include <string>

namespace chimera {

class CHIMERA_IMEX Symbol {
	const std::string	*sym;
public:
		Symbol();
	explicit Symbol(char);
	explicit Symbol(const char *cp);
	explicit Symbol(const std::string &orig);
	bool	operator==(Symbol s) const { return sym == s.sym; }
	bool	operator<(Symbol s) const { return *sym < *s.sym; }
	bool	operator!=(Symbol s) const { return sym != s.sym; }
	bool	operator<=(Symbol s) const { return *sym <= *s.sym; }
	bool	operator>(Symbol s) const { return *sym > *s.sym; }
	bool	operator>=(Symbol s) const { return *sym >= *s.sym; }
	bool	empty() const { return sym->empty(); }
	const std::string	&str() const { return *sym; }
	std::string::size_type	size() const { return sym->size(); }
	const char	*c_str() const { return sym->c_str(); }
	unsigned int	hash() const;

	// Init class allows static Symbol's by guaranting that the
	// symbol table exists before the Symbol's are created
	class CHIMERA_IMEX Init {
		static int count;
	public:
		Init();
		~Init();
	};
};

inline bool operator==(Symbol s, const std::string &str)
{
	return s.str() == str;
}

inline bool operator<(Symbol s, const std::string &str)
{
	return s.str() < str;
}

inline bool operator!=(Symbol s, const std::string &str)
{
	return s.str() != str;
}

inline bool operator<=(Symbol s, const std::string &str)
{
	return s.str() <= str;
}

inline bool operator>(Symbol s, const std::string &str)
{
	return s.str() > str;
}

inline bool operator>=(Symbol s, const std::string &str)
{
	return s.str() >= str;
}

inline bool operator==(const std::string &str, Symbol s)
{
	return str == s.str();
}

inline bool operator<(const std::string &str, Symbol s)
{
	return str < s.str();
}

inline bool operator!=(const std::string &str, Symbol s)
{
	return str != s.str();
}

inline bool operator<=(const std::string &str, Symbol s)
{
	return str <= s.str();
}

inline bool operator>(const std::string &str, Symbol s)
{
	return str > s.str();
}

inline bool operator>=(const std::string &str, Symbol s)
{
	return str >= s.str();
}

inline std::ostream &operator<<(std::ostream &os, Symbol s)
{
	return os << s.str();
}

} // namespace chimera

namespace {
// Guarantee that Symbol table is initialized before used.
// This isn't needed for shared libraries that are loaded after
// the chimera library has been loaded.
chimera::Symbol::Init __chimeraSymbolInit;
}

# ifdef CHIMERA_STD_HASH
#  include <hash_set>

template <> struct std::hash<chimera::Symbol> {
	unsigned int operator()(chimera::Symbol s) const
	{
		return s.hash();
	}
};
# endif

#endif
