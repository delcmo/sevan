#ifndef SEVANAPP_H
#define SEVANAPP_H

#include "MooseApp.h"

class SevanApp;

template<>
InputParameters validParams<SevanApp>();

class SevanApp : public MooseApp
{
public:
  SevanApp(const std::string & name, InputParameters parameters);
  virtual ~SevanApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* SEVANAPP_H */
