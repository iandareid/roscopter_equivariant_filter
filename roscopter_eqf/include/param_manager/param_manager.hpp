#ifndef PARAM_MANAGER_H
#define PARAM_MANAGER_H

#include <map>
#include <rclcpp/rclcpp.hpp>
#include <string>
#include <variant>

namespace roscopter
{

class ParamManager
{
public:
  explicit ParamManager(rclcpp::Node * node);

  double get_double(std::string param_name);
  bool get_bool(std::string param_name);
  int64_t get_int(std::string param_name);
  std::string get_string(std::string param_name);

  void declare_double(std::string param_name, double value);
  void declare_bool(std::string param_name, bool value);
  void declare_int(std::string param_name, int64_t value);
  void declare_string(std::string param_name, std::string value);

  void set_double(std::string param_name, double value);
  void set_bool(std::string param_name, bool value);
  void set_int(std::string param_name, int64_t value);
  void set_string(std::string param_name, std::string value);

  void set_parameters();
  bool set_parameters_callback(const std::vector<rclcpp::Parameter> & parameters);

private:
  using ParameterValue = std::variant<double, bool, int64_t, std::string>;

  template<typename T>
  T get_value(std::string param_name);

  std::map<std::string, ParameterValue> params_;
  rclcpp::Node * container_node_;
};

}  // namespace roscopter

#endif  // PARAM_MANAGER_H
