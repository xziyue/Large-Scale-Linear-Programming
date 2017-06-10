#ifndef DEF_ONDISKMATRIXTYPEINFO_HPP
#define DEF_ONDISKMATRIXTYPEINFO_HPP


#include "__include.hpp"


struct TypeInfoDefined {
	static constexpr bool pass = true;
	static constexpr char* info = "";
};

struct TypeInfoUndefined {
	static constexpr bool pass = false;
	static constexpr char* info = "Type info undefined";
};

template<typename CheckClass>
struct TemplateChecker {
	using CheckClassType = CheckClass;
	TemplateChecker() {
		if (!CheckClass::pass) {
			cerr << CheckClass::info << "\n";
			throw runtime_error{ CheckClass::info };
		}
	}
};



template<typename V, typename InfoDefineChecker = TypeInfoUndefined>
struct TypeInfo :public TemplateChecker<InfoDefineChecker> {
	const int32_t type_size = sizeof(V);
	const char type_hint[3] = { '\0','\0','\0' };
};

template<typename InfoDefineChecker>
struct TypeInfo<double, InfoDefineChecker> :public TemplateChecker<TypeInfoDefined> {
	const int32_t type_size = sizeof(double);
	const char type_hint[3] = { 'f','6','4' };
};

template<typename InfoDefineChecker>
struct TypeInfo<float, InfoDefineChecker> :public TemplateChecker<TypeInfoDefined> {
	const int32_t type_size = sizeof(float);
	const char type_hint[3] = { 'f','3','2' };
};


#endif // !DEF_ONDISKMATRIXTYPEINFO_HPP

